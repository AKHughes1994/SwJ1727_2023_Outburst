#!/usr/bin/env python
"""
Iterate over epoch result directories, select models per epoch, and write
per-epoch summary text files.
"""

import json
import re
import sys
from datetime import datetime
from itertools import product
from pathlib import Path

from astropy.time import Time

# =============================================================================
# Configuration
# =============================================================================

ITERATE_ALL = True   # False → process TEST_DIR only

_REPO_ROOT = Path(__file__).parent.parent

RESULTS_DIR = _REPO_ROOT / 'results' / 'J1727'
QU_TEXT_DIR = _REPO_ROOT / 'files'   / 'J1727' / 'QU_text'
SUMMARY_DIR = _REPO_ROOT / 'files'   / 'J1727' / 'summary'
TEST_DIR = RESULTS_DIR / 'SwiftJ1727_WAPITI_20231006'

# =============================================================================
# Complexity sequence
# =============================================================================

def _build_complexity_sequence(max_components=5):
    """
    Return ordered list of valid non-P model types.
    Rules: at least one S, no pure-T models.
    Within a given component count, more S = less complex.
    """
    sequence = []
    for n in range(1, max_components + 1):
        for n_s in range(n, 0, -1):   # n_s from n down to 1
            n_t = n - n_s
            sequence.append('S' * n_s + 'T' * n_t)
    return sequence


COMPLEXITY_SEQUENCE = _build_complexity_sequence()

# =============================================================================
# Helpers
# =============================================================================

_DATE_RE = re.compile(r'(\d{8})')


def _extract_date_string(name):
    """Return the 8-digit date string from a directory or file name."""
    m = _DATE_RE.search(name)
    return m.group(1) if m else None


def _find_qu_text_file(date_str):
    """Find the QU/WAPITI text file matching an 8-digit date string."""
    for p in QU_TEXT_DIR.iterdir():
        if p.suffix == '.txt' and date_str in p.name:
            return p
    return None


def _parse_isot_from_header(txt_path):
    """Read ISOT timestamp from the file header (line 2)."""
    with open(txt_path) as f:
        for line in f:
            m = re.search(r'Date \(ISO\):\s*(\S+)', line)
            if m:
                return m.group(1)
    raise ValueError(f"No ISO date found in {txt_path}")


def _load_json(path):
    with open(path) as f:
        return json.load(f)


def _parse_model_type(filename):
    """
    Extract model type from filename like:
    SwiftJ1727_QU_20231112_PP_results.json → PP
    """
    stem = Path(filename).stem                  # strip .json
    stem = stem.replace('_results', '')         # strip _results
    parts = stem.split('_')
    return parts[-1]                            # last token is model type


def _is_p_model(model_type):
    return bool(re.fullmatch(r'P+', model_type))


def _is_valid_non_p(model_type):
    """At least one S, no pure T, only S and T characters."""
    if not re.fullmatch(r'[ST]+', model_type):
        return False
    if 'S' not in model_type:
        return False
    if re.fullmatch(r'T+', model_type):
        return False
    return True


def _is_pure_s(model_type):
    return bool(re.fullmatch(r'S+', model_type))


def _best_by_evidence(results_map):
    """Return model_type with highest logz from a dict {model_type: json_data}."""
    return max(results_map, key=lambda k: results_map[k]['evidence']['logz'])


# =============================================================================
# Parameter extraction
# =============================================================================

def _extract_processed_params(data):
    """
    Extract mode and 68% HDI for every entry in processed_parameters.
    Returns list of dicts with keys: name, mode, hdi_lo, hdi_hi.
    """
    out = []
    for name, info in data.get('processed_parameters', {}).items():
        mode_block = info.get('mode_1')
        if mode_block is None:
            continue
        hdi = mode_block.get('hdi_68', [None, None])
        out.append({
            'name': name,
            'mode': mode_block.get('mode'),
            'hdi_lo': hdi[0],
            'hdi_hi': hdi[1],
        })
    return out


def _extract_stokes_i(data):
    si = data.get('stokes_i_fit')
    if si is None:
        return None
    return {
        'model':           si.get('model'),
        'I0':              si.get('I0'),
        'I0_err':          si.get('I0_err'),
        'alpha_lambda2':   si.get('alpha_lambda2'),
        'alpha_lambda2_err': si.get('alpha_lambda2_err'),
        'alpha_nu':        si.get('alpha_nu'),
        'alpha_nu_err':    si.get('alpha_nu_err'),
        'freq0_GHz':       si.get('freq0_GHz'),
    }



    gof = data['goodness_of_fit']
    best = gof['best_sample']
    iw   = gof['importance_weighted']
    return {
        'chi2_best':         best['chi_squared'],
        'chi2_best_Q':       best.get('chi_squared_Q'),
        'chi2_best_U':       best.get('chi_squared_U'),
        'reduced_chi2_best': best['reduced_chi_squared'],
        'chi2_mean':         iw['chi_squared_mean'],
        'chi2_std':          iw['chi_squared_std'],
        'dof':               gof['dof'],
    }


def _extract_evidence(data):
    ev = data['evidence']
    return {'logz': ev['logz'], 'logzerr': ev['logzerr']}


def _extract_gof(data):
    gof = data['goodness_of_fit']
    best = gof['best_sample']
    iw   = gof['importance_weighted']
    return {
        'chi2_best':         best['chi_squared'],
        'chi2_best_Q':       best.get('chi_squared_Q'),
        'chi2_best_U':       best.get('chi_squared_U'),
        'reduced_chi2_best': best['reduced_chi_squared'],
        'chi2_mean':         iw['chi_squared_mean'],
        'chi2_std':          iw['chi_squared_std'],
        'dof':               gof['dof'],
    }


def _extract_metadata(data):
    md = data['metadata']
    return {
        'model_type':  md['model_type'],
        'n_params':    md['n_params'],
        'n_data':      md['n_data'],
        'n_iter':      md.get('n_iter'),
        'n_eff':       md.get('n_eff'),
        'bic':         data['model_comparison']['BIC'],
        'aic':         data['model_comparison']['AIC'],
    }


# =============================================================================
# Output formatting
# =============================================================================

_SEP  = '=' * 72
_SEP2 = '-' * 72


def _fmt_param(p):
    lo  = p['hdi_lo']
    hi  = p['hdi_hi']
    md  = p['mode']
    err_lo = (md - lo) if (md is not None and lo is not None) else None
    err_hi = (hi - md) if (md is not None and hi is not None) else None

    def _f(v):
        return f'{v:+.6f}' if v is not None else 'N/A'

    return (
        f"  {p['name']:<30s}"
        f"mode = {_f(md)}   "
        f"-{abs(err_lo):.6f} / +{abs(err_hi):.6f}  (68% HDI)"
        if err_lo is not None and err_hi is not None
        else f"  {p['name']:<30s}mode = {_f(md)}   HDI unavailable"
    )


def _write_model_block(lines, label, data):
    md  = _extract_metadata(data)
    ev  = _extract_evidence(data)
    gof = _extract_gof(data)
    si  = _extract_stokes_i(data)
    params = _extract_processed_params(data)

    lines.append(_SEP2)
    lines.append(f"  {label}: {md['model_type']}")
    lines.append(_SEP2)
    lines.append(f"  n_params = {md['n_params']}   "
                 f"n_data = {md['n_data']}   "
                 f"DOF = {gof['dof']}   "
                 f"n_iter = {md['n_iter']}   "
                 f"n_eff = {md['n_eff']:.0f}" if md['n_eff'] else '')
    lines.append(f"  BIC = {md['bic']:.3f}   AIC = {md['aic']:.3f}")
    lines.append('')
    lines.append('  Evidence')
    lines.append(f"    ln(Z) = {ev['logz']:.4f} ± {ev['logzerr']:.4f}")
    lines.append('')
    lines.append('  Goodness of fit')
    lines.append(f"    chi2 (best sample) = {gof['chi2_best']:.3f}   "
                 f"[Q={gof['chi2_best_Q']:.3f}  U={gof['chi2_best_U']:.3f}]"
                 if gof['chi2_best_Q'] is not None else
                 f"    chi2 (best sample) = {gof['chi2_best']:.3f}")
    lines.append(f"    reduced chi2 (best) = {gof['reduced_chi2_best']:.4f}")
    lines.append(f"    chi2 (importance-weighted) = "
                 f"{gof['chi2_mean']:.3f} ± {gof['chi2_std']:.3f}")
    lines.append('')
    lines.append('  Parameters  [mode  -lo/+hi  68% HDI]')
    for p in params:
        lines.append(_fmt_param(p))
    lines.append('')
    if si is not None:
        lines.append('  Stokes I fit')
        lines.append(f"    model  = {si['model']}   ν₀ = {si['freq0_GHz']} GHz")
        lines.append(f"    I₀     = {si['I0']:.6f} ± {si['I0_err']:.6f}")
        lines.append(f"    α(λ²)  = {si['alpha_lambda2']:.6f} ± {si['alpha_lambda2_err']:.6f}")
        lines.append(f"    α(ν)   = {si['alpha_nu']:.6f} ± {si['alpha_nu_err']:.6f}")
        lines.append('')


def _write_epoch_summary(epoch_dir, date_isot, date_mjd, selected, out_path):
    """
    selected = {
        'best_s_series':   json_data or None,
        'best_non_p':      json_data,
        'next_complexity': json_data,
        'best_p':          json_data or None,
    }
    """
    lines = []
    lines.append(_SEP)
    lines.append(f"  EPOCH SUMMARY: {epoch_dir.name}")
    lines.append(_SEP)
    lines.append(f"  Date (ISOT) : {date_isot}")
    lines.append(f"  Date (MJD)  : {date_mjd:.6f}")
    lines.append('')

    label_map = [
        ('SINGLE S MODEL',           'single_s'),
        ('BEST S-SERIES MODEL',      'best_s_series'),
        ('BEST NON-P MODEL',         'best_non_p'),
        ('NEXT COMPLEXITY MODEL',    'next_complexity'),
        ('BEST P MODEL',             'best_p'),
    ]
    for label, key in label_map:
        data = selected.get(key)
        if data is None:
            lines.append(_SEP2)
            lines.append(f"  {label}: not available for this epoch")
            lines.append('')
        else:
            _write_model_block(lines, label, data)

    lines.append(_SEP)
    lines.append('')

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text('\n'.join(lines))
    print(f"  Written: {out_path.name}")


# =============================================================================
# Per-epoch processing
# =============================================================================

def process_epoch(epoch_dir):
    date_str = _extract_date_string(epoch_dir.name)
    if date_str is None:
        print(f"  [SKIP] Cannot parse date from: {epoch_dir.name}")
        return

    # Match QU text file and parse date
    txt_file = _find_qu_text_file(date_str)
    if txt_file is None:
        print(f"  [SKIP] No QU text file found for date {date_str}")
        return

    date_isot = _parse_isot_from_header(txt_file)
    date_mjd  = Time(date_isot, format='isot', scale='utc').mjd

    # Load all results JSON files in this directory
    json_files = list(epoch_dir.glob('*_results.json'))
    if not json_files:
        print(f"  [SKIP] No results JSON files in {epoch_dir.name}")
        return

    all_models = {}
    for jf in json_files:
        mtype = _parse_model_type(jf.name)
        all_models[mtype] = _load_json(jf)

    # Partition by category
    p_models      = {k: v for k, v in all_models.items() if _is_p_model(k)}
    pure_s_models = {k: v for k, v in all_models.items() if _is_pure_s(k)}
    non_p_models  = {k: v for k, v in all_models.items() if _is_valid_non_p(k)}

    # Ejecta epochs: single S model only
    if 'ejecta' in epoch_dir.name.lower():
        single_s_data = all_models.get('S', None)
        if single_s_data is None:
            print(f"  [SKIP] Ejecta epoch but no S model found in {epoch_dir.name}")
            return
        selected = {'single_s': single_s_data}
        out_path = SUMMARY_DIR / f"{epoch_dir.name}_summary.txt"
        _write_epoch_summary(epoch_dir, date_isot, date_mjd, selected, out_path)
        return

    # Single S model (always reported if present)
    single_s_data = all_models.get('S', None)

    # Best S-series
    best_s_data = None
    if pure_s_models:
        best_s_key  = _best_by_evidence(pure_s_models)
        best_s_data = pure_s_models[best_s_key]

    # Best non-P (must exist)
    if not non_p_models:
        print(f"  [SKIP] No valid non-P models in {epoch_dir.name}")
        return
    # Build collection of valid non-P models within ΔlnZ ≤ 10 of max,
    # sorted by complexity (position in COMPLEXITY_SEQUENCE)
    max_logz = max(v['evidence']['logz'] for v in non_p_models.values())
    within_threshold = [
        mt for mt in COMPLEXITY_SEQUENCE
        if mt in non_p_models
        and max_logz - non_p_models[mt]['evidence']['logz'] <= 10.0
    ]

    if not within_threshold:
        print(f"  [SKIP] No non-P models within threshold for {epoch_dir.name}")
        return

    best_non_p_key  = within_threshold[0]
    best_non_p_data = non_p_models[best_non_p_key]

    # Next least complex within the same collection (None if only one member)
    next_data = None
    if len(within_threshold) > 1:
        next_key  = within_threshold[1]
        next_data = non_p_models[next_key]

    # Best P: simplest within ΔlnZ ≤ 10 of max
    best_p_data = None
    if p_models:
        max_p_logz = max(v['evidence']['logz'] for v in p_models.values())
        p_sequence = [mt for mt in ['P', 'PP', 'PPP', 'PPPP', 'PPPPP']
                      if mt in p_models
                      and max_p_logz - p_models[mt]['evidence']['logz'] <= 10.0]
        if p_sequence:
            best_p_data = p_models[p_sequence[0]]

    selected = {
        'single_s':        single_s_data,
        'best_s_series':   best_s_data,
        'best_non_p':      best_non_p_data,
        'next_complexity': next_data,
        'best_p':          best_p_data,
    }

    # Remove duplicate model blocks (same model type appearing under multiple labels)
    seen_types = set()
    for key, data in list(selected.items()):
        if data is None:
            continue
        mtype = data['metadata']['model_type']
        if mtype in seen_types:
            selected[key] = None
        else:
            seen_types.add(mtype)

    out_path = SUMMARY_DIR / f"{epoch_dir.name}_summary.txt"
    _write_epoch_summary(epoch_dir, date_isot, date_mjd, selected, out_path)


# =============================================================================
# Main
# =============================================================================

def main():
    if ITERATE_ALL:
        epoch_dirs = sorted([
            d for d in RESULTS_DIR.iterdir()
            if d.is_dir() and 'trials' not in d.name.lower()
        ])
    else:
        epoch_dirs = [TEST_DIR]

    print(f"Processing {len(epoch_dirs)} epoch(s)...\n")

    errors = []
    for epoch_dir in epoch_dirs:
        print(f"[{epoch_dir.name}]")
        try:
            process_epoch(epoch_dir)
        except RuntimeError as e:
            print(f"  [ERROR] {e}")
            errors.append((epoch_dir.name, str(e)))
        except Exception as e:
            print(f"  [ERROR] Unexpected: {e}")
            errors.append((epoch_dir.name, str(e)))

    print(f"\nDone. {len(epoch_dirs) - len(errors)} succeeded, "
          f"{len(errors)} failed.")
    if errors:
        print("\nFailed epochs:")
        for name, msg in errors:
            print(f"  {name}: {msg}")


if __name__ == '__main__':
    main()
