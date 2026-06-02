#!/usr/bin/env python3
"""
xf_batch_export.py
──────────────────
Batch-converts MeerKAT/CASA *.Xf cross-hand phase calibration tables to
plain-text files suitable for manual_XF_solver / oxkat downstream use.

Directory layout assumed
────────────────────────
  XF_ROOT_DIR/
    <any nesting>/
      cal_1GC_<CODE>_<ms_stem>.Xf/        ← CASA table directory
      cal_1GC_<CODE>_<ms_stem>_combineScan.Xf/

Output naming
─────────────
  OUTPUT_DIR/<CODE>_<YYYYMMDD>_<FIELD>_XF.txt

  For tables that contain solutions for more than one field, one file is
  written per field.

Usage
─────
  Adjust XF_ROOT_DIR and OUTPUT_DIR below, then:
      python xf_batch_export.py
  or with a CASA-aware python:
      casa --nologger --nogui -c xf_batch_export.py
"""

from __future__ import annotations

import sys
import re
from pathlib import Path
from datetime import datetime, timezone, timedelta

import numpy as np

try:
    from pyrap.tables import table as casatable
except ImportError:
    # Allow the import to fail at module level so editors can still parse the
    # file; the error will surface clearly at runtime.
    casatable = None  # type: ignore

# ── Global configuration ──────────────────────────────────────────────────────

_HERE = Path(__file__).resolve().parent   # SpecPolJetComp/PlotsAndPostProc/code

# Parent directory that contains (possibly nested) *.Xf CASA table directories.
# NOTE: The raw Xf CASA tables are NOT included in this repository due to GitHub
# storage constraints. The extracted plain-text output files are already provided
# in files/crosshand_phase_files/ — you do not need to re-run this script to use
# them. This script is included only in case it is useful for someone who wants
# to re-extract from the original CASA tables.
XF_ROOT_DIR = _HERE.parent / 'xf_tables'

# Where the output *.txt files will be written (already provided in the repo).
OUTPUT_DIR = _HERE.parent / 'files/crosshand_phase_files'

# Glob pattern used to find Xf tables under XF_ROOT_DIR.
# '**/*.Xf' descends into all sub-directories; change to '*.Xf' for flat layout.
XF_GLOB = '*.Xf'

# ── Filename parsing ──────────────────────────────────────────────────────────

_CODE_RE = re.compile(r'(?:^|_)(\d{10})(?:_|$)')


def extract_code(xfpath: Path) -> str:
    """Return the 10-digit Unix-timestamp code embedded in the Xf filename.

    e.g. cal_1GC_1694011170_sdp_l0_1024ch.ms.Xf  →  '1694011170'
    """
    m = _CODE_RE.search(xfpath.name)
    if m:
        return m.group(1)
    raise ValueError(f"Cannot extract 10-digit code from filename: {xfpath.name}")


def is_combine_scan(xfpath: Path) -> bool:
    return '_combineScan' in xfpath.name


# ── Time conversion ───────────────────────────────────────────────────────────

def mjd_sec_to_utc(mjd_sec_arr: np.ndarray) -> datetime:
    """Convert an array of MJD-seconds (CASA TIME column) to a UTC datetime."""
    mjd_days = float(np.mean(mjd_sec_arr)) / 86400.0
    unix_sec = (mjd_days - 40587.0) * 86400.0   # MJD 40587 = 1970-01-01
    return datetime(1970, 1, 1, tzinfo=timezone.utc) + timedelta(seconds=unix_sec)


# ── Core table reader ─────────────────────────────────────────────────────────

def read_xf_table(xftab: Path):
    """Read all relevant columns from a CASA Xf gain table.

    Returns
    -------
    phase_deg : ndarray, shape (nrows, nchan)
        Phase in degrees for the first polarisation.
    spw_ids   : ndarray, shape (nrows,)
    field_ids : ndarray, shape (nrows,)
    flag      : ndarray or None, shape (nrows, nchan)
    time_mjd  : ndarray, shape (nrows,)   [MJD seconds]
    chan_freqs : ndarray, shape (n_spw, n_chan)
    field_names: list of str
    """
    if casatable is None:
        raise RuntimeError("pyrap.tables is not importable – check your environment.")

    # ── Main table ────────────────────────────────────────────────────────────
    tb       = casatable(str(xftab), readonly=True, ack=False)
    colnames = tb.colnames()

    if 'FPARAM' in colnames:
        param       = tb.getcol('FPARAM')          # (nrows, nchan, npol)  – real radians
        phase_deg   = np.degrees(param[..., 0])
    elif 'CPARAM' in colnames:
        param       = tb.getcol('CPARAM')          # (nrows, nchan, npol)  – complex
        phase_deg   = np.degrees(np.angle(param[..., 0]))
    else:
        tb.close()
        raise RuntimeError(
            f"Neither FPARAM nor CPARAM found in {xftab.name}.\n"
            f"Available columns: {colnames}"
        )

    spw_ids   = tb.getcol('SPECTRAL_WINDOW_ID')    # (nrows,)
    field_ids = tb.getcol('FIELD_ID')               # (nrows,)
    time_mjd  = tb.getcol('TIME')                   # (nrows,)  MJD seconds
    flag      = tb.getcol('FLAG')[..., 0] if 'FLAG' in colnames else None
    tb.close()

    # ── SPW sub-table ─────────────────────────────────────────────────────────
    spw_tb     = casatable(str(xftab) + '/SPECTRAL_WINDOW', readonly=True, ack=False)
    chan_freqs  = spw_tb.getcol('CHAN_FREQ')         # (n_spw, n_chan)
    spw_tb.close()

    # ── FIELD sub-table ───────────────────────────────────────────────────────
    field_names: list[str] = []
    field_path = Path(str(xftab)) / 'FIELD'
    if field_path.is_dir():
        field_tb    = casatable(str(field_path), readonly=True, ack=False)
        raw         = field_tb.getcol('NAME')       # list of strings
        field_names = [str(n).strip() for n in raw]
        field_tb.close()
    else:
        # Fallback: name fields by integer ID
        max_fid = int(field_ids.max()) + 1
        field_names = [f'field{i}' for i in range(max_fid)]

    return phase_deg, spw_ids, field_ids, flag, time_mjd, chan_freqs, field_names


# ── Per-field output writer ───────────────────────────────────────────────────

def process_field(
    xftab: Path,
    field_id: int,
    field_name: str,
    phase_deg: np.ndarray,
    spw_ids: np.ndarray,
    field_ids: np.ndarray,
    flag: np.ndarray | None,
    time_mjd: np.ndarray,
    chan_freqs: np.ndarray,
    obs_dt: datetime,
    code: str,
    output_dir: Path,
) -> Path:
    """Median-collapse time axis for one field and write the output text file."""

    field_mask  = field_ids == field_id
    f_spws      = spw_ids[field_mask]
    f_phases    = phase_deg[field_mask]         # (n_field_rows, n_chan)
    f_flags     = flag[field_mask] if flag is not None else None

    freq_list: list[np.ndarray]  = []
    phase_list: list[np.ndarray] = []

    for spw in np.unique(f_spws):
        rows   = f_spws == spw
        phases = f_phases[rows]                  # (n_times, n_chan)

        if f_flags is not None:
            phases = np.where(f_flags[rows], np.nan, phases)

        med_phase = np.nanmedian(phases, axis=0)  # (n_chan,)
        freqs     = chan_freqs[spw]

        valid = np.isfinite(med_phase)
        freq_list.append(freqs[valid])
        phase_list.append(med_phase[valid])

    freq_out  = np.concatenate(freq_list)
    phase_out = np.concatenate(phase_list)
    order     = np.argsort(freq_out)
    freq_out  = freq_out[order]
    phase_out = phase_out[order]

    date_str   = obs_dt.strftime('%Y%m%d')
    time_str   = obs_dt.strftime('%H:%M:%S')
    out_name   = f"{code}_{date_str}_{field_name}_XF.txt"
    out_path   = output_dir / out_name
    output_dir.mkdir(parents=True, exist_ok=True)

    header = "\n".join([
        "Cross-hand phase (Xf) calibration table export",
        f"# ──────────────────────────────────────────────",
        f"# Source table  : {xftab.name}",
        f"# Code          : {code}",
        f"# Obs date (UTC): {obs_dt.strftime('%Y-%m-%d')}",
        f"# Obs time (UTC): {time_str}",
        f"# Field ID      : {field_id}",
        f"# Field name    : {field_name}",
        f"# CombineScan   : {'yes' if is_combine_scan(xftab) else 'no'}",
        f"# Generated     : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"# N channels    : {len(freq_out)}",
        f"# Freq range    : {freq_out.min()/1e9:.4f} – {freq_out.max()/1e9:.4f} GHz",
        f"# Median phase  : {np.nanmedian(phase_out):.4f} deg",
        f"# Std phase     : {np.nanstd(phase_out):.4f} deg",
        f"#",
        f"# Columns: freq(Hz)   phase(deg)",
    ])

    np.savetxt(out_path, np.column_stack([freq_out, phase_out]),
               header=header, fmt='%.10e', comments='# ')
    return out_path


# ── Main batch loop ───────────────────────────────────────────────────────────

def main() -> None:
    xf_tables = sorted(XF_ROOT_DIR.glob(XF_GLOB))[:]

    if not xf_tables:
        print(f"[WARNING] No *.Xf tables found under {XF_ROOT_DIR}", file=sys.stderr)
        return

    print(f"Found {len(xf_tables)} Xf table(s) under {XF_ROOT_DIR}\n")

    n_ok = n_err = 0

    for xftab in xf_tables:
        print(f"Processing : {xftab.name}")
        try:
            code = extract_code(xftab)

            (phase_deg, spw_ids, field_ids, flag,
             time_mjd, chan_freqs, field_names) = read_xf_table(xftab)

            obs_dt = mjd_sec_to_utc(time_mjd)

            unique_fields = np.unique(field_ids)
            for fid in unique_fields:
                fname = (field_names[fid]
                         if fid < len(field_names)
                         else f'field{fid}')
                # Sanitise field name for use in a filename
                fname_safe = re.sub(r'[^\w\-+.]', '_', fname)

                out_path = process_field(
                    xftab, fid, fname_safe,
                    phase_deg, spw_ids, field_ids, flag,
                    time_mjd, chan_freqs, obs_dt, code, OUTPUT_DIR,
                )
                n_ok += 1
                print(f"  ✓ field '{fname}' → {out_path.name}")

        except Exception as exc:
            n_err += 1
            print(f"  ✗ FAILED: {exc}", file=sys.stderr)

        print()

    print(f"Done. {n_ok} file(s) written, {n_err} error(s).")


if __name__ == '__main__':
    main()
