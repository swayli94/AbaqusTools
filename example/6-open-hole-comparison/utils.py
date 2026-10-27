'''
Utility functions for the comparison of FEA and IM results.
'''
import os
path = os.path.dirname(os.path.abspath(__file__))

import numpy as np
from typing import Dict, Tuple, List
from scipy.interpolate import griddata

N_CASE = 3
N_RADIAL = 32
N_CIRCUMFERENTIAL = 64

NAME_COMPONENTS_PLOT = ['S11', 'S22', 'S12']

SOURCES       = ['fem_C3D8R', 'fem_S4R', 'im_S4R']
SOURCE_LABELS = ['FEM C3D8R', 'FEM S4R', 'IM S4R']
CASE_LABELS   = [f'Case {i}' for i in range(N_CASE)]
N_SRC  = len(SOURCES)
N_COMP = len(NAME_COMPONENTS_PLOT)


def create_unit_template_mesh(r_hole: float = 1.0,
                r_outer: float = 2.0) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Create a unit template mesh for the 2D field reconstruction.
    The mesh is a circular mesh with a hole in the middle (hole radius = 1.0, outer radius = 2.0).
    '''
    r = np.linspace(0, 1, N_RADIAL, endpoint=True)
    r = r_hole + (r_outer - r_hole) * r**2 # Dense the points near the hole boundary.
    theta = np.linspace(0, 2*np.pi, N_CIRCUMFERENTIAL, endpoint=True)
    R, Theta = np.meshgrid(r, theta, indexing='ij')
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)
    return X, Y

def interpolate_to_template(raw_data: np.ndarray,
        X_template: np.ndarray, Y_template: np.ndarray) -> np.ndarray:
    '''
    Interpolate scattered field data (hole-centered) onto the structured template mesh.

    Parameters
    ----------
    raw_data : np.ndarray, shape (n_points, 2 + len(NAME_COMPONENTS))
        Hole-centered coordinates and component values.
    X_template, Y_template : np.ndarray, shape (N_RADIAL, N_CIRCUMFERENTIAL)
        Template mesh coordinates.

    Returns
    -------
    field : np.ndarray, shape (N_RADIAL, N_CIRCUMFERENTIAL, len(NAME_COMPONENTS))
    '''
    xy_src = raw_data[:, :2]
    values_src = raw_data[:, 2:]
    n_comp = values_src.shape[1]
    xy_query = np.column_stack([X_template.ravel(), Y_template.ravel()])

    field = np.zeros((*X_template.shape, n_comp))
    for i in range(n_comp):
        vals = griddata(xy_src, values_src[:, i], xy_query, method='linear')
        # Fill any NaN (outside source convex hull) with nearest-neighbour.
        nan_mask = np.isnan(vals)
        if np.any(nan_mask):
            vals[nan_mask] = griddata(xy_src, values_src[:, i],
                                      xy_query[nan_mask], method='nearest')
        field[:, :, i] = vals.reshape(X_template.shape)

    return field

def load_field_data(z_planes: List[float], i_sample: int, i_case: int,
                source: str = 'fem_C3D8R',
                path_data: str = os.path.join(path, 'data', 'fem_C3D8R')
                ) -> Dict[float, np.ndarray]:
    '''
    Load and extract the field data from each ply z-plane.

    File formats
    ------------
    fem_C3D8R: Job_OHP_{s}_{c}-field.dat
        Variables= X Y Z index S11 S22 S33 S12 S13 S23

    fem_S4R: Job_OHP_{s}_{c}-field.dat
        Variables= X Y Z index S11 S22 S12 TSHR13 TSHR23 E11 E22 E12

    im_S4R: Job_OHP_{s}_{c}-stress-field.dat
        Variables= X Y Z S11 S22 S12

    Returns
    -------
    field_data : Dict[float, np.ndarray]
        A dictionary mapping z-planes to their corresponding field data.
        Each value is an array of shape (n_points, 6).
        Columns: [x, y, z, S11, S22, S12].
        Coordinates are still in the original FEA frame (not hole-centered)
        for fem_C3D8R and fem_S4R; already hole-centered for im_S4R.
    '''
    if source == 'fem_C3D8R':
        # Variables= X Y Z index S11 S22 S33 S12 S13 S23
        fname = f'Job_OHP_{i_sample}_{i_case}-field.dat'
        index_components = [0, 1, 2, 4, 5, 7] # x, y, z, S11, S22, S12
    elif source == 'fem_S4R':
        # Variables= X Y Z index S11 S22 S12 TSHR13 TSHR23 E11 E22 E12
        fname = f'Job_OHP_{i_sample}_{i_case}-field.dat'
        index_components = [0, 1, 2, 4, 5, 6] # x, y, z, S11, S22, S12
    elif source == 'im_S4R':
        # Variables= X Y Z S11 S22 S12
        fname = f'Job_OHP_{i_sample}_{i_case}-stress-field.dat'
        index_components = [0, 1, 2, 3, 4, 5] # x, y, z, S11, S22, S12
    else:
        raise ValueError(f'Unknown source: {source}')

    fname = os.path.join(path_data, fname)
    data = np.loadtxt(fname, skiprows=2)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    data = data[:, index_components] # Keep only the needed components.

    # Extract the data of each z-plane
    z_unique = np.unique(data[:, 2])
    field_data = {}

    for z in z_planes:
        if not np.any(np.isclose(z_unique, z)):
            print(f'>>> Warning: z-plane {z} not found in data for sample {i_sample}, '
                  f'source {source}. Available planes: {z_unique}')
            continue
        field_data[z] = data[np.isclose(data[:, 2], z)]

    return field_data
