'''
Slope-ratio analysis comparing FEM S4R and IM S4R against FEM C3D8R (reference).

For each prediction source x load case x stress component (aggregated over all plies),
fit y = k·x through the origin (OLS), where x = C3D8R reference, y = prediction.

k is the scale-ratio slope:
    k > 1  →  systematic overshoot  (conservative)
    k < 1  →  systematic undershoot (non-conservative)
    k = 1  →  perfect proportionality

sigma_k is the OLS standard error of k (68% confidence interval: k ± sigma_k).

Metrics saved per (source x case x component):
    slope_all        : k fitted on all valid points
    slope_all_sig    : sigma_k for slope_all
    top5pct_slope    : k fitted on top-5% highest-|sigma_ref| elements
    top5pct_slope_sig: sigma_k for top5pct_slope

Output files (./error-analysis/):
    {label}_error_metrics.csv    — slope metrics table
    {label}_scatter_{src}.png    — pred vs ref scatter with fitted lines annotated
    {label}_slope_errorbars.png  — k ± sigma_k error-bar charts grouped by source
'''

import os
import sys
import csv
path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(path)

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from utils import (
    N_CASE, NAME_COMPONENTS_PLOT,
    SOURCES, SOURCE_LABELS, CASE_LABELS, N_SRC, N_COMP,
    create_unit_template_mesh,
)
from compare_ply_stress_field import derive_ply_info, collect_fields

path_out = os.path.join(path, 'error-analysis')
os.makedirs(path_out, exist_ok=True)

DPI   = 150
TOP_K = 5.0   # percentage of top-stress elements used for dangerous-zone slope

def load_rp11_rf(i_sample: int, i_case: int, source: str, path_data: str) -> float:
    '''
    Parse Job_OHP_{i_sample}_{i_case}-RF.dat and return RP_11_RF.
    Returns nan if the file is missing or the key is absent.
    '''
    fname = os.path.join(path_data, f'Job_OHP_{i_sample}_{i_case}-RF.dat')
    if not os.path.exists(fname):
        print(f'  Warning: {fname} not found')
        return np.nan
    with open(fname) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2 and parts[0] == 'RP_11_RF':
                try:
                    return float(parts[1])
                except ValueError:
                    pass
    return np.nan


def collect_rp11_rf(i_sample_dict: dict, func_path_data) -> np.ndarray:
    '''
    Load RP_11_RF for all (case x source) combinations.

    Returns
    -------
    np.ndarray, shape (N_CASE, N_SRC)
    '''
    out = np.full((N_CASE, N_SRC), np.nan)
    for i_case in range(N_CASE):
        for i_src, src in enumerate(SOURCES):
            out[i_case, i_src] = load_rp11_rf(
                i_sample=i_sample_dict[src],
                i_case=i_case,
                source=src,
                path_data=func_path_data(src),
            )
    return out


def plot_rf_comparison(rp11_rf: np.ndarray, label: str):
    '''
    Grouped bar chart: RP_11_RF across sources and cases.
    '''
    x     = np.arange(N_CASE)
    width = 0.8 / N_SRC

    fig, ax = plt.subplots(figsize=(N_CASE * 1.8 + 1.5, 3.5))

    for i_src, src_label in enumerate(SOURCE_LABELS):
        vals   = rp11_rf[:, i_src]
        offset = (i_src - (N_SRC - 1) / 2.0) * width
        ax.bar(x + offset, vals, width=width * 0.92, label=src_label, zorder=2)

    ax.axhline(0.0, color='dimgray', ls='--', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(CASE_LABELS, fontsize=8)
    ax.set_ylabel('RP_11_RF (N)', fontsize=8)
    ax.grid(axis='y', lw=0.4, alpha=0.5, zorder=0)
    ax.legend(fontsize=7)

    fig.suptitle(f'Reaction force RP_11_RF  |  {label}', fontsize=10)
    plt.tight_layout()
    fname = os.path.join(path_out, f'{label}_rf_comparison.png')
    fig.savefig(fname, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {fname}')


def fit_slope_through_origin(x: np.ndarray, y: np.ndarray) -> tuple:
    '''
    OLS fit of y = k*x (no intercept).

    k       = Σ(x_i y_i) / Σ(x_i²)
    sigma_k = sqrt(s² / Σ(x_i²)),  s² = Σ(y_i - k x_i)² / (n-1)

    Returns (k, sigma_k); both nan when data are insufficient.
    k=1 → perfect proportionality; k>1 → overshoot; k<1 → undershoot.
    '''
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = x.size
    if n < 2:
        return np.nan, np.nan
    xx = float(np.dot(x, x))
    if xx < 1e-30:
        return np.nan, np.nan
    k  = float(np.dot(x, y) / xx)
    s2 = float(np.dot(y - k * x, y - k * x)) / max(n - 1, 1)
    return k, float(np.sqrt(s2 / xx))


def compute_sample_metrics(fields: np.ndarray, k: float = TOP_K) -> list:
    '''
    Compute slope metrics for every (pred_source x case x component) triple.

    fields : shape (N_CASE, N_SRC, n_ply, NR, NC, N_COMP)
             source index 0 is the C3D8R reference.
    Returns a list of dicts, one per combination.
    '''
    records  = []
    col_topk = f'top{int(k)}pct'

    for i_src in range(1, N_SRC):
        for i_case in range(N_CASE):
            for i_comp, comp in enumerate(NAME_COMPONENTS_PLOT):
                ref  = fields[i_case, 0,     :, :, :, i_comp].ravel()
                pred = fields[i_case, i_src, :, :, :, i_comp].ravel()
                mask = np.isfinite(ref) & np.isfinite(pred)
                r, p = ref[mask], pred[mask]

                thr = np.percentile(np.abs(r), 100.0 - k) if r.size else 0.0
                hi  = np.abs(r) >= thr

                slope_all, slope_all_sig = fit_slope_through_origin(r,     p)
                slope_top, slope_top_sig = fit_slope_through_origin(r[hi], p[hi])

                records.append({
                    'source':                SOURCE_LABELS[i_src],
                    'case':                  CASE_LABELS[i_case],
                    'component':             comp,
                    'slope_all':             slope_all,
                    'slope_all_sig':         slope_all_sig,
                    f'{col_topk}_slope':     slope_top,
                    f'{col_topk}_slope_sig': slope_top_sig,
                })
    return records


def save_metrics_csv(records: list, label: str) -> str:
    fname = os.path.join(path_out, f'{label}_error_metrics.csv')
    if not records:
        return fname
    keys = list(records[0].keys())
    with open(fname, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for rec in records:
            row = {k: (f'{v:.5f}' if isinstance(v, float) else v)
                   for k, v in rec.items()}
            w.writerow(row)
    return fname


def plot_scatter_topk(fields: np.ndarray, label: str, k: float = TOP_K):
    '''
    One figure per prediction source.
    Layout: N_CASE rows x N_COMP columns.
    Grey points = rest of the field; red points = top-k% dangerous elements.
    Black dashed line = ideal (y=x); blue dash-dot = slope fitted on all points;
    red dash-dot = slope fitted on top-k% points.
    '''
    for i_src in range(1, N_SRC):
        fig, axes = plt.subplots(N_CASE, N_COMP,
                                 figsize=(N_COMP * 3.2, N_CASE * 3.2),
                                 squeeze=False)

        for i_case in range(N_CASE):
            for i_comp, comp in enumerate(NAME_COMPONENTS_PLOT):
                ax   = axes[i_case][i_comp]
                ref  = fields[i_case, 0,     :, :, :, i_comp].ravel()
                pred = fields[i_case, i_src, :, :, :, i_comp].ravel()
                mask = np.isfinite(ref) & np.isfinite(pred)
                ref, pred = ref[mask], pred[mask]

                thr = np.percentile(np.abs(ref), 100.0 - k)
                hi  = np.abs(ref) >= thr

                ax.scatter(ref[~hi], pred[~hi], s=1.5, alpha=0.25,
                           c='lightgrey', zorder=1, rasterized=True)
                ax.scatter(ref[hi],  pred[hi],  s=4,   alpha=0.7,
                           c='tab:red', zorder=2, label=f'top {k:.0f}%',
                           rasterized=True)

                lv = max(float(np.percentile(np.abs(np.r_[ref, pred]), 99)), 1e-6)
                ax.plot([-lv, lv], [-lv, lv], 'k--', lw=0.8, label='ideal (k=1)')

                xr = np.array([-lv, lv])
                k_all, sk_all = fit_slope_through_origin(ref,     pred)
                k_top, sk_top = fit_slope_through_origin(ref[hi], pred[hi])
                if np.isfinite(k_all):
                    ax.plot(xr, k_all * xr, color='steelblue',
                            lw=1.0, ls='-.', alpha=0.85, zorder=3)
                if np.isfinite(k_top):
                    ax.plot(xr, k_top * xr, color='tab:red',
                            lw=1.2, ls='-.', alpha=0.85, zorder=3)

                ann_lines = []
                if np.isfinite(k_all):
                    ann_lines.append(f'k_all   = {k_all:.2f} ± {sk_all:.3f}')
                if np.isfinite(k_top):
                    ann_lines.append(f'k_top{int(k):.0f}% = {k_top:.2f} ± {sk_top:.3f}')
                if ann_lines:
                    ax.text(0.03, 0.97, '\n'.join(ann_lines),
                            transform=ax.transAxes, fontsize=10.0,
                            va='top', ha='left', family='monospace',
                            bbox=dict(boxstyle='round,pad=0.25',
                                      fc='white', alpha=0.82, lw=0))

                ax.set_xlim(-lv * 1.12, lv * 1.12)
                ax.set_ylim(-lv * 1.12, lv * 1.12)
                ax.set_xlabel(f'{comp} C3D8R (MPa)', fontsize=8)
                ax.set_ylabel(f'{comp} {SOURCE_LABELS[i_src]} (MPa)', fontsize=8)
                ax.set_title(f'{CASE_LABELS[i_case]} / {comp}', fontsize=8)
                if i_case == 0 and i_comp == N_COMP - 1:
                    ax.legend(fontsize=7, markerscale=2)

        fig.suptitle(f'{SOURCE_LABELS[i_src]} vs FEM C3D8R  |  {label}', fontsize=10)
        plt.tight_layout()
        fname = os.path.join(path_out, f'{label}_scatter_{SOURCES[i_src]}.png')
        fig.savefig(fname, dpi=DPI, bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved {fname}')


def plot_slope_errorbars(records: list, label: str, k: float = TOP_K):
    '''
    Single figure with 2 rows x N_COMP columns.
    Row 0: k fitted on all points.
    Row 1: k fitted on top-k% points.
    Each subplot shows k ± sigma_k as error bars, grouped by source, over cases.
    Horizontal dashed line at k=1 (ideal).
    '''
    col_topk  = f'top{int(k)}pct'
    row_specs = [
        ('slope_all',         'slope_all_sig',         'All points'),
        (f'{col_topk}_slope', f'{col_topk}_slope_sig', f'Top-{int(k)}%'),
    ]

    sources = list(dict.fromkeys(r['source'] for r in records))
    cases   = list(dict.fromkeys(r['case']   for r in records))
    x       = np.arange(len(cases))
    width   = 0.8 / len(sources)

    n_rows = len(row_specs)
    fig, axes = plt.subplots(n_rows, N_COMP,
                             figsize=(N_COMP * 4.0, n_rows * 3.2),
                             squeeze=False)

    for i_row, (col_mean, col_sig, row_title) in enumerate(row_specs):
        for i_comp, comp in enumerate(NAME_COMPONENTS_PLOT):
            ax = axes[i_row][i_comp]

            for j, src in enumerate(sources):
                means, sigs = [], []
                for c in cases:
                    hit = next((r for r in records
                                if r['source'] == src and r['case'] == c
                                and r['component'] == comp), None)
                    means.append(hit[col_mean] if hit else np.nan)
                    sigs.append(hit[col_sig]   if hit else np.nan)
                offset = (j - (len(sources) - 1) / 2.0) * width
                ax.errorbar(x + offset, means, yerr=sigs,
                            fmt='o', capsize=4, capthick=1.2,
                            markersize=5, linewidth=1.2, label=src)

            ax.axhline(1.0, color='dimgray', ls='--', lw=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(cases, fontsize=8)
            ax.set_ylabel('k  (ideal = 1)', fontsize=8)
            ax.set_title(f'{row_title} — {comp}', fontsize=8)
            if i_row == 0 and i_comp == N_COMP - 1:
                ax.legend(fontsize=7)

    fig.suptitle(f'Slope-ratio k ± σ  |  {label}', fontsize=10)
    plt.tight_layout()
    fname = os.path.join(path_out, f'{label}_slope_errorbars.png')
    fig.savefig(fname, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {fname}')


if __name__ == '__main__':

    path_example = os.path.dirname(path)
    path_root = os.path.dirname(path_example)

    PATH_DATA = {
        'fem_C3D8R': os.path.join(path_example, '6-open-hole-specimen-C3D8R', 'data'),
        'fem_S4R':   os.path.join(path_example, '6-open-hole-specimen-S4R',   'data'),
        'im_S4R':    os.path.join(path_example, '6-open-hole-implicit-modelling-S4R', 'data'),
        # 'im_S4R':    os.path.join(path_root, 'temp-MIF-S4R'),
    }

    fname_params = os.path.join(PATH_DATA['fem_C3D8R'], '..', 'default-parameters.json')
    with open(fname_params) as f:
        params = json.load(f)
    pGeo  = params['pGeo']
    pMesh = params['pMesh']

    i_sample = 0
    n_ply, ply_orientations, z_planes = derive_ply_info(pMesh)

    r_hole = pGeo['r_hole']
    X_unit, Y_unit = create_unit_template_mesh(r_hole=1.0, r_outer=1.5)
    X_tmpl = X_unit * r_hole
    Y_tmpl = Y_unit * r_hole

    print(f'Sample {i_sample}: r_hole={r_hole} mm, '
          f'n_ply={n_ply}, orientations={ply_orientations}')

    i_sample_dict = {src: i_sample for src in SOURCES}

    fields = collect_fields(i_sample_dict, pGeo, n_ply, z_planes, X_tmpl, Y_tmpl,
                            func_path_data=lambda src: PATH_DATA[src])

    sample_label = f'sample{i_sample}'

    records  = compute_sample_metrics(fields)
    csv_path = save_metrics_csv(records, sample_label)
    print(f'  Saved {csv_path}')

    plot_scatter_topk(fields, sample_label)
    plot_slope_errorbars(records, sample_label)

    rp11_rf = collect_rp11_rf(i_sample_dict,
                              func_path_data=lambda src: PATH_DATA[src])
    plot_rf_comparison(rp11_rf, sample_label)
