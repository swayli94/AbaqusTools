'''
Compare the ply-level stress field of different modelling approaches:
- explicit modelling of bolts with C3D8R elements (7-bolted-joint-specimen-C3D8R)
- explicit modelling of bolts with SC8R elements (7-bolted-joint-specimen-SC8R)
- mesh-independent fastener modelling with S4R elements (7-mesh-independent-fastener-S4R)

The calculated results are stored in the 'data' folder.

Comparison of stress components: S11, S22, S12
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
from scipy.interpolate import griddata
from typing import Callable, Dict, List, Tuple

N_CASE            = 3
N_RADIAL          = 32
N_CIRCUMFERENTIAL = 64
R_OUTER_RATIO     = 1.4    # outer template radius in units of r_hole

NAME_COMPONENTS_PLOT = ['S11', 'S22', 'S12']

SOURCES       = ['fem_C3D8R', 'fem_SC8R', 'mif_S4R',         'mif_S4R_IM']
SOURCE_LABELS = ['FEM C3D8R', 'FEM SC8R', 'MIF S4R (no IM)', 'MIF S4R (IM)']
CASE_LABELS   = [f'Case {i}' for i in range(N_CASE)]
N_SRC  = len(SOURCES)
N_COMP = len(NAME_COMPONENTS_PLOT)

DPI = 150

# Mesh utilities

def create_unit_template_mesh(r_hole: float = 1.0,
                              r_outer: float = R_OUTER_RATIO
                              ) -> Tuple[np.ndarray, np.ndarray]:
    '''Annular template mesh; r is denser near the hole boundary.'''
    r = np.linspace(0, 1, N_RADIAL, endpoint=True)
    r = r_hole + (r_outer - r_hole) * r**2
    theta = np.linspace(0, 2 * np.pi, N_CIRCUMFERENTIAL, endpoint=True)
    R, Theta = np.meshgrid(r, theta, indexing='ij')
    return R * np.cos(Theta), R * np.sin(Theta)

def interpolate_to_template(raw_data: np.ndarray,
                             X_template: np.ndarray,
                             Y_template: np.ndarray) -> np.ndarray:
    '''
    Interpolate scattered (hole-centred) data onto the structured template mesh.

    raw_data : (n_pts, 2 + n_comp)  — first two cols are x, y
    Returns field of shape (N_RADIAL, N_CIRCUMFERENTIAL, n_comp).
    '''
    xy_src     = raw_data[:, :2]
    values_src = raw_data[:, 2:]
    n_comp     = values_src.shape[1]
    xy_query   = np.column_stack([X_template.ravel(), Y_template.ravel()])

    field = np.zeros((*X_template.shape, n_comp))
    for i in range(n_comp):
        vals = griddata(xy_src, values_src[:, i], xy_query, method='linear')
        nan_mask = np.isnan(vals)
        if np.any(nan_mask):
            vals[nan_mask] = griddata(xy_src, values_src[:, i],
                                      xy_query[nan_mask], method='nearest')
        field[:, :, i] = vals.reshape(X_template.shape)
    return field

# Data loading

def _parse_dat_rows(fname: str, n_col: int) -> list:
    '''Read a Tecplot .dat file, returning only rows with exactly n_col floats.'''
    rows = []
    with open(fname) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('Variables') or s.startswith('zone'):
                continue
            try:
                vals = [float(x) for x in s.split()]
                if len(vals) == n_col:
                    rows.append(vals)
            except ValueError:
                pass
    return rows

def load_field_data(z_planes: List[float], i_sample: int, i_case: int,
                    source: str, path_data: str,
                    plate_idx: int = 0,
                    plate_thickness: float = 5.0) -> Dict[float, np.ndarray]:
    '''
    Load ply-level stress data for one (sample, case, source, plate).

    Parameters
    ----------
    plate_idx       : 0 = PLATE_0 (bottom), 1 = PLATE_1 (top)
    plate_thickness : total thickness of one plate (mm), used to split FEM data

    File formats
    ------------
    fem_C3D8R: Job_BJ_{s}_{c}-field.dat
        Variables= X Y Z index S11 S22 S33 S12 S13 S23
        PLATE_0: Z ∈ [0, plate_thickness); PLATE_1: Z ∈ [plate_thickness, 2*plate_thickness)

    fem_SC8R: Job_BJ_{s}_{c}-field-SC8R.dat
        Variables= X Y Z thickness S11 S22 S33 S12 index index_thickness
        Z (col 2) is the element centroid — used to identify the plate.
        "thickness" (col 3) is the ply mid-plane z used as the z-plane identifier.

    mif_S4R, mif_S4R_IM: Job_MIF_{s}_{c}-stress-field-{plate_idx}.dat
        Variables= X Y Z S11 S22 S12
        Coordinates are already hole-centred; one file per plate.

    Returns
    -------
    Dict mapping ply z-plane (mm) → array of shape (n_pts, 6):
        [x, y, z_ply_local, S11, S22, S12]
    z_ply_local is normalized to [0, plate_thickness) for all sources.
    For fem_C3D8R / fem_SC8R, x and y are in the global FEA frame.
    For mif_S4R and mif_S4R_IM, x and y are already hole-centred.
    '''
    z_lo = plate_idx * plate_thickness
    z_hi = (plate_idx + 1) * plate_thickness

    if source == 'fem_C3D8R':
        # Variables= X Y Z index S11 S22 S33 S12 S13 S23
        fname = os.path.join(path_data, f'Job_BJ_{i_sample}_{i_case}-field.dat')
        rows  = _parse_dat_rows(fname, n_col=10)
        arr   = np.array(rows)
        arr   = arr[(arr[:, 2] >= z_lo) & (arr[:, 2] < z_hi)]  # keep plate
        arr   = arr[:, [0, 1, 2, 4, 5, 7]]                     # X Y Z S11 S22 S12
        arr[:, 2] -= z_lo                                       # normalize to local z

    elif source == 'fem_SC8R':
        # Variables= X Y Z thickness S11 S22 S33 S12 index index_thickness
        fname = os.path.join(path_data, f'Job_BJ_{i_sample}_{i_case}-field-SC8R.dat')
        rows  = _parse_dat_rows(fname, n_col=10)
        arr   = np.array(rows)
        arr   = arr[(arr[:, 2] >= z_lo) & (arr[:, 2] < z_hi)]  # filter by element Z
        arr   = arr[:, [0, 1, 3, 4, 5, 7]]                     # X Y thickness S11 S22 S12
        # thickness is already local to each plate — no z offset needed

    elif source in ['mif_S4R', 'mif_S4R_IM']:
        # Variables= X Y Z S11 S22 S12  (one file per plate, already hole-centred)
        fname = os.path.join(path_data,
                             f'Job_MIF_{i_sample}_{i_case}-stress-field-{plate_idx}.dat')
        if not os.path.exists(fname):
            raise FileNotFoundError(
                f'No mif_S4R stress-field file for sample {i_sample} '
                f'case {i_case} plate {plate_idx}')
        rows = _parse_dat_rows(fname, n_col=6)
        arr  = np.array(rows)    # X Y Z S11 S22 S12

    else:
        raise ValueError(f'Unknown source: {source}')

    # Column 2 is the ply z-coordinate for all three sources after remapping above.
    z_col      = arr[:, 2]
    field_data = {}
    for z in z_planes:
        mask = np.isclose(z_col, z, atol=1e-3)
        if not np.any(mask):
            print(f'  >>> Warning: z={z:.3f} mm not found '
                  f'[{source} s{i_sample} c{i_case}]')
            continue
        field_data[z] = arr[mask]
    return field_data

# Ply geometry

def derive_ply_info(pMesh: dict) -> tuple:
    '''
    Compute ply count, orientations, and mid-plane z-coordinates from pMesh.

    Returns
    -------
    n_ply : int
    ply_orientations : list[int]   full stack, bottom → top
    z_planes : list[float]         mid-plane z of each ply (mm)
    '''
    half = pMesh['plate_CompositePly_orientationValue']
    full = list(half) + list(reversed(half)) \
           if pMesh['plate_CompositeLayup_symmetric'] else list(half)
    t    = pMesh['composite_ply_thickness']
    n    = len(full)
    return n, full, [(i + 0.5) * t for i in range(n)]

# Field collection

def collect_fields(i_sample_dict: Dict[str, int],
                   pGeo: dict, n_ply: int, z_planes: list,
                   X_tmpl: np.ndarray, Y_tmpl: np.ndarray,
                   func_path_data: Callable[[str], str],
                   plate_idx: int = 0,
                   plate_thickness: float = 5.0) -> np.ndarray:
    '''
    Load and interpolate stress fields for all cases, sources, and plies.

    Parameters
    ----------
    plate_idx       : which plate to load (0 or 1)
    plate_thickness : total thickness of one plate (mm)

    Returns
    -------
    np.ndarray, shape (N_CASE, N_SRC, n_ply, N_RADIAL, N_CIRCUMFERENTIAL, N_COMP)
    '''
    fastener = pGeo['fasteners'][0]
    cx = fastener['x_center']
    cy = fastener['y_center']

    out = np.full((N_CASE, N_SRC, n_ply, N_RADIAL, N_CIRCUMFERENTIAL, N_COMP), np.nan)

    for i_case in range(N_CASE):
        for i_src, src in enumerate(SOURCES):
            try:
                fd = load_field_data(z_planes,
                                     i_sample=i_sample_dict[src],
                                     i_case=i_case, source=src,
                                     path_data=func_path_data(src),
                                     plate_idx=plate_idx,
                                     plate_thickness=plate_thickness)
            except Exception as e:
                print(f'  Warning [{src} s{i_sample_dict[src]} c{i_case}]: {e}')
                continue

            for i_ply, z in enumerate(z_planes):
                if z not in fd:
                    continue
                raw = fd[z]          # (n, 6): x  y  z_ply  S11  S22  S12
                xy  = raw[:, :2].copy()
                if src not in ['mif_S4R', 'mif_S4R_IM']:  # FEM data is in the global plate frame
                    xy[:, 0] -= cx
                    xy[:, 1] -= cy
                out[i_case, i_src, i_ply] = interpolate_to_template(
                    np.column_stack([xy, raw[:, 3:6]]), X_tmpl, Y_tmpl
                )

    return out

# Plot functions

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
    Layout: (N_CASE x N_SRC) rows x n_ply columns.
    All rows in the same case block share one symmetric colorbar.
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
                    _draw_field(ax, X, Y,
                                fields[i_case, i_src, i_ply, :, :, i_comp],
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

        fig.suptitle(f'{comp} Stress Field — {sample_label}', fontsize=10)
        fname = os.path.join(path_out, f'{sample_label}_{comp}_field.png')
        fig.savefig(fname, dpi=DPI, bbox_inches='tight')
        plt.close(fig)


if __name__ == '__main__':

    path_example = os.path.dirname(path)
    path_root = os.path.dirname(path_example)

    PATH_DATA = {
        'fem_C3D8R':  os.path.join(path_example, '7-bolted-joint-specimen-C3D8R', 'data'),
        'fem_SC8R':   os.path.join(path_example, '7-bolted-joint-specimen-SC8R',  'data'),
        'mif_S4R':    os.path.join(path_example, '7-mesh-independent-fastener-S4R', 'data', 'no-IM'),
        'mif_S4R_IM': os.path.join(path_example, '7-mesh-independent-fastener-S4R', 'data', 'with-IM'),
        # 'mif_S4R_IM': os.path.join(path_root, 'temp-MIF-S4R'),
    }

    path_figure = os.path.join(path, 'field-figure')
    os.makedirs(path_figure, exist_ok=True)

    fname_params = os.path.join(PATH_DATA['fem_C3D8R'], '..', 'default-parameters.json')
    with open(fname_params) as f:
        params = json.load(f)
    pGeo  = params['pGeo']
    pMesh = params['pMesh']

    i_sample = 0
    n_ply, ply_orientations, z_planes = derive_ply_info(pMesh)
    plate_thickness = n_ply * pMesh['composite_ply_thickness']

    r_hole = pGeo['fasteners'][0]['r_hole']
    X_unit, Y_unit = create_unit_template_mesh(r_hole=1.0, r_outer=R_OUTER_RATIO)
    X_tmpl = X_unit * r_hole
    Y_tmpl = Y_unit * r_hole

    print(f'Sample {i_sample}: r_hole={r_hole} mm, '
          f'n_ply={n_ply}, '
          f'plate_thickness={plate_thickness} mm')

    i_sample_dict = {src: i_sample for src in SOURCES}

    for plate_idx in range(2):
        fields = collect_fields(i_sample_dict, pGeo, n_ply, z_planes, X_tmpl, Y_tmpl,
                                func_path_data=lambda src: PATH_DATA[src],
                                plate_idx=plate_idx,
                                plate_thickness=plate_thickness)
        sample_label = f'sample{i_sample}_plate{plate_idx}'
        plot_component_figures(fields, X_tmpl, Y_tmpl, path_figure, sample_label,
                               ply_orientations)
