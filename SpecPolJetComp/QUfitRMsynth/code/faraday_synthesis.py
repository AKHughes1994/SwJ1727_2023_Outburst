import warnings
warnings.filterwarnings('ignore', message='.*pkg_resources.*')

import re
import numpy as np
from pathlib import Path
from collections import defaultdict
try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        raise ImportError("Python ≥3.11 required for tomllib, or install tomli: pip install tomli")
from faraday_plot import plot_flagging_diagnostic
from faraday_data import DataLoader, FlaggingConfig
import subprocess

# Flagging parameters shared with faraday_obs_fit.py; single source of truth
FLAGGING_CONFIG_FILE = Path(__file__).parent / 'flagging_config.toml'

# Targets to process; each must have a corresponding files/<TARGET>/QU_text/ directory
TARGETS = ['J1727', '3C286', 'J1733']

# Per-epoch overrides for rmclean thresholds.
# These epochs have bright, complex FDFs that need a deeper clean (lower peak threshold)
# and a tighter finer-mask to avoid cleaning into noise.
# Keyed by frozenset of date substrings so multiple epochs share one entry.
PEAK_OVERRIDES = {
    frozenset(['1014', '1006', '1016']): '-14',
}
FINER_OVERRIDES = {
    frozenset(['1014', '1006', '1016']): '-6',
}

# Fallback thresholds used when no override matches.
# -9 ≈ 10^-9 × peak → very deep clean; -4 ≈ 10^-4 × peak for the finer mask.
DEFAULT_PEAK  = '-9'
DEFAULT_FINER = '-4'


def get_peak_threshold(base_name):
    """Return the rmclean1d peak threshold (-c) for this epoch."""
    for keys, val in PEAK_OVERRIDES.items():
        if any(k in base_name for k in keys):
            return val
    return DEFAULT_PEAK


def get_finer_threshold(base_name):
    """Return the rmclean1d finer-mask threshold (-w) for this epoch."""
    for keys, val in FINER_OVERRIDES.items():
        if any(k in base_name for k in keys):
            return val
    return DEFAULT_FINER


def _date_key(stem):
    """Extract YYYYMMDD from a file stem for grouping; fall back to the full stem."""
    m = re.search(r'(\d{8})', stem)
    return m.group(1) if m else stem


def _resolve_files(files_path):
    """
    Return the ordered list of QU text files to process.

    J1733 has multiple scans per epoch stored as separate files (e.g. scan08, scan10).
    For any date where scan files exist, process each scan individually so that
    per-scan RM synthesis can later be averaged or inspected.
    For all other targets/dates, process one file per epoch as usual.
    """
    all_files = sorted(files_path.glob('*.txt'))
    groups = defaultdict(list)
    for f in all_files:
        groups[_date_key(f.stem)].append(f)

    resolved = []
    for date in sorted(groups):
        group = groups[date]
        scan_files = [f for f in group if re.search(r'scan', f.stem, re.IGNORECASE)]
        if len(group) > 1 and scan_files:
            print(f"  [{date}] {len(scan_files)} scan file(s) found — running per-scan")
            resolved.extend(sorted(scan_files))
        else:
            resolved.extend(group)
    return resolved


def process_target(target, flagging_cfg):
    files_path = Path(__file__).parent.parent / 'files' / target / 'QU_text'
    output_path = Path(__file__).parent.parent / 'rmsynth' / target
    output_path.mkdir(parents=True, exist_ok=True)   # create rmsynth/<target>/ if absent
    input_files = _resolve_files(files_path)

    if not input_files:
        print(f"No input files found for {target}, skipping.")
        return

    # Check once up-front whether any epoch has already been processed, and if so
    # ask the user whether to skip or redo those epochs.  Avoids re-running hours
    # of RM synthesis unless explicitly requested.
    any_existing = any(
        (output_path / f.stem / f"{f.stem}_rmsynth_FDFmodel.dat").exists()
        for f in input_files
    )
    refit_existing = False
    if any_existing:
        response = input(f"\nSome outputs already exist for {target}.\nRe-fit already processed files? (y/n): ").strip().lower()
        refit_existing = response == 'y'

    for input_file in input_files:
        print(f"\nProcessing: {input_file}")
        base_name = input_file.stem
        output_subdir = output_path / base_name
        output_subdir.mkdir(parents=True, exist_ok=True)  # one subdirectory per epoch/scan

        # Skip if already processed and user did not ask to refit
        already_fit = (output_subdir / f"{base_name}_rmsynth_FDFmodel.dat").exists()
        if already_fit and not refit_existing:
            print(f"Skipping {base_name} (already fit).")
            continue

        # Stokes I median is used to set the jitter floor in the flagging config —
        # channels whose residual from the spectral fit exceeds jitter_frac × median(I)
        # are flagged as bad.
        stokes_i = np.median(np.loadtxt(input_file)[:, 1])

        # Apply the default flagging config, then overwrite with the first matching
        # per-epoch override (if any).  The override system mirrors faraday_obs_fit.py
        # so that the same channels are masked in both the RM synthesis and the QU fitting.
        fparams = dict(flagging_cfg['default'])
        for override in flagging_cfg.get('overrides', []):
            if any(m in base_name for m in override['match']):
                fparams.update({k: v for k, v in override.items() if k != 'match'})
                print(f"Flagging: applied override (match={override['match']})")
                break
        else:
            print("Flagging: using default configuration")

        spectrum_jitter = fparams['jitter_frac'] * stokes_i
        print(f"Stokes I median: {stokes_i:.6f}, jitter: {spectrum_jitter:.6f}")

        flagging = FlaggingConfig(
            sigma_clip_extreme=fparams['sigma_clip_extreme'],
            freq_ranges=fparams['freq_ranges'],
            clip_stokes_i=None,
            max_iter=100,
            sigma_clip_errors=fparams['sigma_clip_errors'],
            sigma_clip_errors_space='log',
            sigma_clip_spectrum=fparams['sigma_clip_spectrum'],
            spectrum_deg=fparams['spectrum_deg'],
            spectrum_jitter=spectrum_jitter,
            spectrum_clip_sided=fparams['spectrum_clip_sided'],
            spectrum_stokes=fparams['spectrum_stokes'],
            sigma_clip_v=fparams['sigma_clip_v'],
            clip_v_against='median',
            verbose=True,
        )

        load_result = DataLoader.load_physical_with_fit(input_file, I_model='polynomial',
                                                        flagging=flagging, load_stokes_v=True)
        data, I_fitter, flagging_info = load_result

        # Save a diagnostic plot showing which channels were flagged and why
        plot_flagging_diagnostic(flagging_info['raw_data'], flagging_info,
                                 filename=f'{base_name}_flagging_diagnostic.png',
                                 output_dir=output_subdir)

        # Write the flagged physical-unit spectrum to disk.
        # rmsynth1d expects columns: freq[Hz]  I  Q  U  dI  dQ  dU
        # We use the raw (pre-fractional) values from flagging_info so that
        # rmsynth1d can apply its own noise weighting.
        flagged_cols = [flagging_info['raw_data'][k][~flagging_info['mask']]
                        for k in ['freq', 'I', 'Q', 'U', 'dI', 'dQ', 'dU']]
        rmsynth_file = output_subdir / f"{base_name}_rmsynth.txt"
        np.savetxt(rmsynth_file, np.column_stack(flagged_cols))

        # Run RM synthesis via CIRADA RM-Tools.
        # -S: save intermediate products
        # -o 4: oversampling factor (4× the nominal channel width)
        # -l 2500: maximum Faraday depth in rad/m²
        # --super-resolution: use super-resolution RM synthesis for finer φ sampling
        subprocess.run([
            'rmsynth1d',
            str(rmsynth_file),
            '-S', '-o', '4', '-l', '2500', '--super-resolution'
        ])

        # Run RM-CLEAN on the dirty FDF produced above.
        # -S: save output
        # -c: log10 of the peak fractional threshold at which to stop cleaning
        # -w: log10 of the finer-mask fractional threshold (tighter mask for final iterations)
        # Thresholds are set per-epoch via PEAK_OVERRIDES / FINER_OVERRIDES above.
        subprocess.run([
            'rmclean1d',
            str(rmsynth_file),
            '-S', '-c', get_peak_threshold(base_name), '-w', get_finer_threshold(base_name)
        ])

def main():
    with open(FLAGGING_CONFIG_FILE, 'rb') as f:
        flagging_cfg = tomllib.load(f)

    for target in TARGETS:
        print(f"\n{'='*60}\nTarget: {target}\n{'='*60}")
        process_target(target, flagging_cfg)


if __name__ == "__main__":
    main()
