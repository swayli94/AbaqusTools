'''
Compare the ply-level stress field of different modelling approaches:
- explicit FEM with C3D8R elements (6-open-hole-specimen-C3D8R)
- explicit FEM with S4R elements  (6-open-hole-specimen-S4R)
- implicit modelling with S4R elements (6-open-hole-implicit-modelling-S4R)

The calculated results are stored in the 'data' folder of each specimen directory.

Comparison of stress components: S11, S22, S12

Three figure types are produced per sample:
  {label}_{comp}_field.png   — ply-by-ply stress field, grouped by case
  {label}_envelope.png       — max(|stress|) envelope across all plies
  {label}_{comp}_ratio.png   — log-ratio (source / C3D8R) field, grouped by case
'''

import os
import sys
path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(path)

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from typing import Callable, Dict

from utils import (
    N_CASE, N_RADIAL, N_CIRCUMFERENTIAL, NAME_COMPONENTS_PLOT,
    SOURCES, SOURCE_LABELS, CASE_LABELS, N_SRC, N_COMP,
    load_field_data,
    create_unit_template_mesh, interpolate_to_template,
)

path_figure = os.path.join(path, 'field-figure')
os.makedirs(path_figure, exist_ok=True)


def derive_ply_info(pMesh: dict) -> tuple:
    '''
    Compute ply count, orientations, and mid-plane z-coordinates from pMesh.

    Returns
    -------
    n_ply : int
    ply_orientations : list[int]   full stack (bottom → top)
    z_planes : list[float]         mid-plane z of each ply (mm)
    '''
    half = pMesh['plate_CompositePly_orientationValue']
    if pMesh['plate_CompositeLayup_symmetric']:
        full = list(half) + list(reversed(half))
    else:
        full = list(half)
    t = pMesh['composite_ply_thickness']
    n = len(full)
    z_planes = [(i + 0.5) * t for i in range(n)]
    return n, full, z_planes


def collect_fields(i_sample_dict: Dict[str, int],
                   pGeo: dict, n_ply: int, z_planes: list,
                   X_tmpl: np.ndarray, Y_tmpl: np.ndarray,
                   func_path_data: Callable[[str], str]) -> np.ndarray:
    '''
    Load and interpolate stress fields for all cases, sources, and plies.

    Returns
    -------
    np.ndarray of shape (N_CASE, N_SRC, n_ply, N_RADIAL, N_CIRCUMFERENTIAL, N_COMP)
    '''
    cx = pGeo['len_x_plate'] * pGeo['xr_hole_center']
    cy = pGeo['len_y_plate'] * pGeo['yr_hole_center']
    out = np.full((N_CASE, N_SRC, n_ply, N_RADIAL, N_CIRCUMFERENTIAL, N_COMP), np.nan)

    for i_case in range(N_CASE):
        for i_src, src in enumerate(SOURCES):
            try:
                fd = load_field_data(z_planes,
                        i_sample=i_sample_dict[src],
                        i_case=i_case, source=src,
                        path_data=func_path_data(src))
            except Exception as e:
                print(f'  Warning [{src} s{i_sample_dict[src]} c{i_case}]: {e}')
                continue

            for i_ply, z in enumerate(z_planes):
                if z not in fd:
                    continue
                raw = fd[z]          # (n, 6): x  y  z  S11  S22  S12
                xy  = raw[:, :2].copy()
                if src != 'im_S4R':  # FEM data is in the global plate frame
                    xy[:, 0] -= cx
                    xy[:, 1] -= cy
                field = interpolate_to_template(
                    np.column_stack([xy, raw[:, 3:6]]), X_tmpl, Y_tmpl
                )
                out[i_case, i_src, i_ply] = field

    return out


def _draw_field(ax, X, Y, data, norm, cmap):
    ax.pcolormesh(X, Y, data, cmap=cmap, norm=norm, shading='gouraud', rasterized=True)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def _add_colorbar(fig, ax_group, norm, cmap, label):
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax_group, fraction=0.025, pad=0.03, shrink=0.85)
    cb.set_label(label, fontsize=7)
    cb.ax.tick_params(labelsize=6)


def plot_component_figures(fields: np.ndarray, X: np.ndarray, Y: np.ndarray,
                           path_out: str, sample_label: str,
                           ply_orientations: list):
    '''
    One figure per stress component.
    Layout: 9 rows (3 cases x 3 sources) x n_ply columns.
    Subplots belonging to the same case share the same colorbar.
    '''
    cmap   = 'RdBu_r'
    n_rows = N_CASE * N_SRC
    n_ply  = fields.shape[2]

    for i_comp, comp in enumerate(NAME_COMPONENTS_PLOT):
        fig, axes = plt.subplots(n_rows, n_ply,
                                 figsize=(n_ply * 1.6 + 2.0, n_rows * 1.6 + 0.8),
                                 squeeze=False)

        for i_case in range(N_CASE):
            block = fields[i_case, :, :, :, :, i_comp]
            valid = block[np.isfinite(block)]
            vlim  = float(np.percentile(np.abs(valid), 99)) if valid.size else 1.0
            vlim  = max(vlim, 1e-6)
            norm  = mcolors.Normalize(vmin=-vlim, vmax=vlim)

            for i_src in range(N_SRC):
                row = i_case * N_SRC + i_src
                for i_ply in range(n_ply):
                    ax = axes[row, i_ply]
                    _draw_field(ax, X, Y, fields[i_case, i_src, i_ply, :, :, i_comp],
                                norm, cmap)
                    if i_ply == 0:
                        ax.text(-0.08, 0.5,
                                f'{CASE_LABELS[i_case]}\n{SOURCE_LABELS[i_src]}',
                                transform=ax.transAxes, fontsize=6,
                                ha='right', va='center')
                    if row == 0:
                        ax.set_title(f'P{i_ply}\n{ply_orientations[i_ply]}°', fontsize=6)

            _add_colorbar(fig,
                          axes[i_case * N_SRC:(i_case + 1) * N_SRC, :],
                          norm, cmap, f'{comp} (MPa)')

        fig.suptitle(f'{comp} Stress Field — sample {sample_label}', fontsize=10)
        fname = os.path.join(path_out, f'{sample_label}_{comp}_field.png')
        fig.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved {fname}')


def plot_envelope_figure(fields: np.ndarray, X: np.ndarray, Y: np.ndarray,
                         path_out: str, sample_label: str):
    '''
    Envelope figure.
    Layout: 9 rows (3 cases x 3 sources) x 3 data columns (components).
    Each cell shows max(|stress|) across all plies.
    Each (case, component) block of 3 rows shares one colorbar placed in a
    dedicated narrow column next to the data column.
    '''
    cmap   = 'viridis'
    n_rows = N_CASE * N_SRC

    # GridSpec: interleave data columns (width 3) and colorbar columns (width 0.15)
    width_ratios = []
    for _ in range(N_COMP):
        width_ratios.extend([3.0, 0.15])

    fig = plt.figure(figsize=(N_COMP * 3.5 + 1.5, n_rows * 2.0 + 0.8))
    gs  = fig.add_gridspec(n_rows, N_COMP * 2,
                           width_ratios=width_ratios,
                           wspace=0.05, hspace=0.08,
                           left=0.10, right=0.97, top=0.95, bottom=0.02)

    # Pre-create data axes (even columns)
    data_axes = [[fig.add_subplot(gs[row, i_comp * 2])
                  for i_comp in range(N_COMP)]
                 for row in range(n_rows)]

    for i_case in range(N_CASE):
        for i_comp, comp in enumerate(NAME_COMPONENTS_PLOT):
            envs = [
                np.nanmax(np.abs(fields[i_case, i_src, :, :, :, i_comp]), axis=0)
                for i_src in range(N_SRC)
            ]
            all_v = np.concatenate([e[np.isfinite(e)] for e in envs])
            vlim  = float(np.percentile(all_v, 99)) if all_v.size else 1.0
            vlim  = max(vlim, 1e-6)
            norm  = mcolors.Normalize(vmin=0, vmax=vlim)

            for i_src, env in enumerate(envs):
                row = i_case * N_SRC + i_src
                ax  = data_axes[row][i_comp]
                _draw_field(ax, X, Y, env, norm, cmap)
                if i_comp == 0:
                    ax.text(-0.12, 0.5,
                            f'{CASE_LABELS[i_case]}\n{SOURCE_LABELS[i_src]}',
                            transform=ax.transAxes, fontsize=7,
                            ha='right', va='center')
                if row == 0:
                    ax.set_title(f'|{comp}|\nenvelope', fontsize=8)

            # One colorbar per (case, component) block, spanning 3 rows
            cbar_ax = fig.add_subplot(
                gs[i_case * N_SRC:(i_case + 1) * N_SRC, i_comp * 2 + 1]
            )
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cb = fig.colorbar(sm, cax=cbar_ax)
            cb.set_label(f'|{comp}| (MPa)', fontsize=7)
            cb.ax.tick_params(labelsize=6)

    fig.suptitle(f'Stress Envelope — sample {sample_label}', fontsize=10)
    fname = os.path.join(path_out, f'{sample_label}_envelope.png')
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {fname}')


def plot_ratio_figures(fields: np.ndarray, X: np.ndarray, Y: np.ndarray,
                       path_out: str, sample_label: str,
                       ply_orientations: list):
    '''
    One figure per stress component.
    Layout: 9 rows (3 cases x 3 sources) x n_ply data columns + 1 colorbar column.
    - fem_C3D8R row : signed original stress field.
    - fem_S4R / im_S4R rows: log-ratio ln(source/C3D8R), masked where |C3D8R| is small.
    Each case has separate colorbars for the reference row and the ratio rows.
    '''
    THRESHOLD = 0.2   # mask ratio where |C3D8R| < 20% of its per-case max
    cmap      = 'RdBu_r'
    n_rows    = N_CASE * N_SRC
    n_ply     = fields.shape[2]

    # Build display array.
    # C3D8R rows keep the original signed field.
    # S4R / IM rows: log-ratio = ln(source / C3D8R), NaN where |C3D8R| is below threshold.
    display = fields.copy()
    for i_case in range(N_CASE):
        for i_comp in range(N_COMP):
            ref      = fields[i_case, 0, :, :, :, i_comp]   # (n_ply, NR, NC)
            ref_safe = ref.copy()
            ref_safe[np.abs(ref) < THRESHOLD * np.nanmax(np.abs(ref))] = np.nan
            for i_src in range(1, N_SRC):
                ratio = fields[i_case, i_src, :, :, :, i_comp] / ref_safe
                log_ratio = np.where(ratio > 0, np.log(ratio), np.nan)
                log_ratio[np.isnan(log_ratio)] = 0.0  # neutral: ln(1)=0
                display[i_case, i_src, :, :, :, i_comp] = log_ratio

    ROW_LABELS = [
        SOURCE_LABELS[0],
        f'{SOURCE_LABELS[1]}\n/ C3D8R',
        f'{SOURCE_LABELS[2]}\n/ C3D8R',
    ]

    for i_comp, comp in enumerate(NAME_COMPONENTS_PLOT):
        norms_ref   = []
        norms_ratio = []
        for i_case in range(N_CASE):
            ref_data = fields[i_case, 0, :, :, :, i_comp]
            valid    = ref_data[np.isfinite(ref_data)]
            vlim     = float(np.percentile(np.abs(valid), 99)) if valid.size else 1.0
            vlim     = max(vlim, 1e-6)
            norms_ref.append(mcolors.Normalize(vmin=-vlim, vmax=vlim))

            ratio_data = display[i_case, 1:, :, :, :, i_comp]
            valid_r    = ratio_data[np.isfinite(ratio_data)]
            if valid_r.size:
                vmax = max(abs(float(np.percentile(valid_r, 1))),
                           abs(float(np.percentile(valid_r, 99))))
            else:
                vmax = np.log(2)
            vmax = max(vmax, 0.01)
            norms_ratio.append(mcolors.Normalize(vmin=-vmax, vmax=vmax))

        width_ratios = [1.0] * n_ply + [0.5]
        fig = plt.figure(figsize=(n_ply * 1.6 + 2.5, n_rows * 1.6 + 0.8))
        gs  = fig.add_gridspec(n_rows, n_ply + 1,
                               width_ratios=width_ratios,
                               wspace=0.05, hspace=0.08,
                               left=0.10, right=0.97,
                               top=0.95, bottom=0.02)

        data_axes = [[fig.add_subplot(gs[row, i_ply])
                      for i_ply in range(n_ply)]
                     for row in range(n_rows)]

        cbar_ref_axes = [
            fig.add_subplot(gs[i_case * N_SRC, n_ply])
            for i_case in range(N_CASE)
        ]
        cbar_ratio_axes = [
            fig.add_subplot(gs[i_case * N_SRC + 1:(i_case + 1) * N_SRC, n_ply])
            for i_case in range(N_CASE)
        ]

        for i_case in range(N_CASE):
            for i_src in range(N_SRC):
                row  = i_case * N_SRC + i_src
                norm = norms_ref[i_case] if i_src == 0 else norms_ratio[i_case]
                for i_ply in range(n_ply):
                    ax = data_axes[row][i_ply]
                    _draw_field(ax, X, Y,
                                display[i_case, i_src, i_ply, :, :, i_comp],
                                norm, cmap)
                    if i_ply == 0:
                        ax.text(-0.08, 0.5,
                                f'{CASE_LABELS[i_case]}\n{ROW_LABELS[i_src]}',
                                transform=ax.transAxes, fontsize=6,
                                ha='right', va='center')
                    if row == 0:
                        ax.set_title(f'P{i_ply}\n{ply_orientations[i_ply]}°', fontsize=6)

        for cbar_ax, norm in zip(cbar_ref_axes, norms_ref):
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cb = fig.colorbar(sm, cax=cbar_ax)
            cb.set_label(f'{comp} (MPa)', fontsize=7)
            cb.ax.tick_params(labelsize=6)

        for cbar_ax, norm in zip(cbar_ratio_axes, norms_ratio):
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cb = fig.colorbar(sm, cax=cbar_ax)
            cb.set_label('ln(source/C3D8R)', fontsize=7)
            cb.ax.tick_params(labelsize=6)

        fig.suptitle(f'{comp} Ratio Field — sample {sample_label}', fontsize=10)
        fname = os.path.join(path_out, f'{sample_label}_{comp}_ratio.png')
        fig.savefig(fname, dpi=150, bbox_inches='tight')
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
    plot_component_figures(fields, X_tmpl, Y_tmpl, path_figure, sample_label,
                           ply_orientations)
    # plot_envelope_figure(fields, X_tmpl, Y_tmpl, path_figure, sample_label)
    # plot_ratio_figures(fields, X_tmpl, Y_tmpl, path_figure, sample_label,
    #                    ply_orientations)
