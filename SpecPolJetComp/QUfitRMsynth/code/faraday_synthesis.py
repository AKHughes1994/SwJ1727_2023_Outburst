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

FLAGGING_CONFIG_FILE = Path(__file__).parent / 'flagging_config.toml'

TARGETS = ['J1727', '3C286', 'J1733']

# rmclean peak and finer-mask thresholds; first match wins
PEAK_OVERRIDES = {
    frozenset(['1014', '1006', '1016']): '-14',
}
FINER_OVERRIDES = {
    frozenset(['1014', '1006', '1016']): '-6',
}

DEFAULT_PEAK  = '-9'
DEFAULT_FINER = '-4'


def get_peak_threshold(base_name):
    for keys, val in PEAK_OVERRIDES.items():
        if any(k in base_name for k in keys):
            return val
    return DEFAULT_PEAK


def get_finer_threshold(base_name):
    for keys, val in FINER_OVERRIDES.items():
        if any(k in base_name for k in keys):
            return val
    return DEFAULT_FINER


def _date_key(stem):
    """Extract YYYYMMDD from a file stem; fall back to the full stem."""
    m = re.search(r'(\d{8})', stem)
    return m.group(1) if m else stem


def _resolve_files(files_path):
    """
    Return the ordered list of files to process.

    For dates where multiple files exist and at least one has 'scan'
    (case-insensitive) in the stem, keep every scan file individually.
    For all other dates, keep the file(s) as-is.
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
    input_files = _resolve_files(files_path)

    if not input_files:
        print(f"No input files found for {target}, skipping.")
        return

    # Ask once if any already-fit outputs exist
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
        output_subdir.mkdir(parents=True, exist_ok=True)

        already_fit = (output_subdir / f"{base_name}_rmsynth_FDFmodel.dat").exists()
        if already_fit and not refit_existing:
            print(f"Skipping {base_name} (already fit).")
            continue

        stokes_i = np.median(np.loadtxt(input_file)[:, 1])

        # Build flagging params: default + first matching override
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
        plot_flagging_diagnostic(flagging_info['raw_data'], flagging_info,
                                 filename=f'{base_name}_flagging_diagnostic.png',
                                 output_dir=output_subdir)

        # Save flagged data for RMSynthesis
        flagged_cols = [flagging_info['raw_data'][k][~flagging_info['mask']]
                        for k in ['freq', 'I', 'Q', 'U', 'dI', 'dQ', 'dU']]
        rmsynth_file = output_subdir / f"{base_name}_rmsynth.txt"
        np.savetxt(rmsynth_file, np.column_stack(flagged_cols))

        # RMSynthesis
        subprocess.run([
            'rmsynth1d',
            str(rmsynth_file),
            '-S', '-o', '4', '-l', '2500', '--super-resolution'
        ])

        # RMClean
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
