'''
Reorganize the field data.

Source data (Tecplot format):
- 3D solid element (C3D8R): X Y Z S11 S22 S33 S12 S13 S23 index
- Continuum shell element (SC8R): X Y Z thickness S11 S22 S12 index index_thickness
- Classical shell element (S4R): X Y Z thickness S11 S22 S12 index index_thickness

The shell element data has three points per ply, the solid element data only has the central point of each ply (element center).
We only need the central point of each ply for all data sources.
For shell elements, we set Z=thickness.

First, we need to match the elements between the data sources,
because the elements with the same coordinates have different indexes (position in the file).

Next, we treat the data points with the same x,y-coordinates as a group,
corresponding to the thickness-distributed ply data of the laminate (shell) element.
Also sort the data points in the ascending order of z-coordinates within each group.

Then, we construct a Pandas DataFrame with the following columns:
- X: x-coordinate of the element center (float)
- Y: y-coordinate of the element center (float)
- total_thickness: laminate thickness (float)
- ply_thickness: thickness of each ply (float)
- n_ply: number of plies (int)
- index: index of the element (int)
- Z-list: list of z-coordinates of the ply centers (list[float])
- S11-list: list of S11 values (list[float])
- S22-list: list of S22 values (list[float])
- S33-list: list of S33 values (list[float])
- S12-list: list of S12 values (list[float])
- S13-list: list of S13 values (list[float])
- S23-list: list of S23 values (list[float])

'''
import os
import sys

path = os.path.dirname(os.path.abspath(__file__))
# Add project root to Python path for multi-branch development
project_root = os.path.abspath(os.path.join(path, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import pandas as pd
import numpy as np
import random

N_COL = 9

fname_C3D8R = 'specimen-stress-field-C3D8R.dat'
fname_SC8R = 'specimen-stress-field-SC8R.dat'
fname_S4R = 'specimen-stress-field-S4R.dat'

fname_C3D8R = os.path.join(path, 'data', fname_C3D8R)
fname_SC8R = os.path.join(path, 'data', fname_SC8R)
fname_S4R = os.path.join(path, 'data', fname_S4R)


def read_tecplot_data(fname: str) -> np.ndarray:
    '''
    Read the Tecplot data file.
    For the shell element data, we only keep the central point of each ply.
    
    Parameters
    ----------
    fname: str
        the file name of the Tecplot data file
    
    Returns
    -------
    data: ndarray [n_point, N_COL]
        the data points: 'X', 'Y', 'Z', 'S11', 'S22', 'S33', 'S12', 'S13', 'S23'
    '''
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    # Find the variables line
    var_line = None
    data_start_idx = None
    for i, line in enumerate(lines):
        if line.startswith('Variables='):
            var_line = line.strip()
            # Find the zone line
            for j in range(i+1, len(lines)):
                if lines[j].strip().startswith('zone'):
                    data_start_idx = j + 1
                    break
            break
    
    if var_line is None or data_start_idx is None:
        raise ValueError(f"Invalid Tecplot format in {fname}")
    
    # Parse variables
    vars_str = var_line.replace('Variables=', '').strip()
    variables = [v.strip() for v in vars_str.split()]
    
    # Determine data format
    is_shell = 'thickness' in variables

    # Read data
    data_lines = []
    for line in lines[data_start_idx:]:
        line = line.strip()
        if not line:
            continue
        # Split by whitespace
        values = line.split()
        if len(values) == 0:
            continue
        data_lines.append(values)
    
    n_point = len(data_lines)
    if n_point == 0:
        raise ValueError(f"No data found in {fname}")

    # Assemble the data points
    if not is_shell:
        # C3D8R source data: X Y Z S11 S22 S33 S12 S13 S23 index
        data = np.array([[float(v) for v in row[:N_COL]] for row in data_lines])
        
    else:
        # SC8R/S4R source data: X Y Z thickness S11 S22 S12 index index_thickness
        data = np.zeros((n_point//3, N_COL))
        i_data = 0
        for i_point in range(n_point):
            
            # Only keep the central point of each ply
            # i.e., index_thickness%3 == 1
            if int(data_lines[i_point][8]) % 3 != 1:
                continue
            
            data[i_data, 0] = float(data_lines[i_point][0]) # X
            data[i_data, 1] = float(data_lines[i_point][1]) # Y
            data[i_data, 2] = float(data_lines[i_point][3]) # Z=thickness
            data[i_data, 3] = float(data_lines[i_point][4]) # S11
            data[i_data, 4] = float(data_lines[i_point][5]) # S22
            data[i_data, 6] = float(data_lines[i_point][6]) # S12
            i_data += 1
    
    return data


def sort_data_points(data: np.ndarray) -> np.ndarray:
    '''
    Sort the data points in the ascending order of x,y,z-coordinates.
    
    Parameters
    ----------
    data: ndarray [n_point, N_COL]
        the data points: 'X', 'Y', 'Z', 'S11', 'S22', 'S33', 'S12', 'S13', 'S23'
    
    Returns
    -------
    data: ndarray [n_point, N_COL]
        the sorted data points
    '''
    # Sort by X, then Y, then Z
    if len(data) == 0:
        return data
    
    # Get indices that would sort the array
    sort_indices = np.lexsort((data[:, 2], data[:, 1], data[:, 0]))
    return data[sort_indices]


def save_as_data_frame(data: np.ndarray, fname_save: str|None = None,
                        tolerance: float = 1e-2) -> pd.DataFrame:
    '''
    Save the data as a Pandas DataFrame.
    If fname_save is not None, save the DataFrame to a CSV file.
    
    Parameters
    ----------
    data: ndarray [n_point, N_COL]
        the data points: 'X', 'Y', 'Z', 'S11', 'S22', 'S33', 'S12', 'S13', 'S23'
    fname_save: str | None, optional
        the file name to save the DataFrame
    tolerance: float, optional
        the tolerance for grouping the data points
        
    Returns
    -------
    df: pd.DataFrame
        the DataFrame of the data points
    '''
    if len(data) == 0:
        return pd.DataFrame()
    
    n_decimal_places = int(np.ceil(np.log10(1/tolerance)))
    
    # Group data by (X, Y) coordinates using quantized coordinates for O(n) grouping
    grouped_data = {}
    
    # Quantize coordinates: round to tolerance precision for fast dictionary lookup
    # This converts coordinates like (20.14672, 39.12858) to quantized keys
    # Points within tolerance will map to the same quantized key
    # Use integer keys to avoid floating point precision issues
    def quantize_coord_to_int(coord: float, tol: float) -> int:
        return int(round(coord / tol))
    
    quantized_to_group = {}
    
    for i, row in enumerate(data):
        x, y = row[0], row[1]
        # Quantize coordinates to integers
        quantized_x_int = quantize_coord_to_int(x, tolerance)
        quantized_y_int = quantize_coord_to_int(y, tolerance)
        quantized_int = (quantized_x_int, quantized_y_int)
        
        if quantized_int in quantized_to_group:
            grouped_data[quantized_int].append(row)
        else:
            quantized_x = round(quantized_x_int * tolerance, n_decimal_places)
            quantized_y = round(quantized_y_int * tolerance, n_decimal_places)
            quantized_to_group[quantized_int] = (quantized_x, quantized_y)
            grouped_data[quantized_int] = [row]
    
    # Check if the number of data points in each group is the same
    n_data_point_in_group = 0
    for quantized_int, points in grouped_data.items():
        if n_data_point_in_group == 0:
            n_data_point_in_group = len(points)
        elif len(points) != n_data_point_in_group:
            print(f'n_data_point_in_group: {n_data_point_in_group}')
            print(f'len(points) in group {quantized_int}: {len(points)}')
            raise ValueError(f"The number of data points in group {quantized_int} is not the same")
    
    if n_data_point_in_group == 0:
        raise ValueError(f"No data points in any group")
    else:
        print(f'  Number of data points in each group, i.e., number of plies: {n_data_point_in_group}')
    
    # Build DataFrame
    # Sort groups by X, then Y coordinates for consistent ordering
    sorted_groups = sorted(grouped_data.items(), key=lambda item: (quantized_to_group[item[0]][0], quantized_to_group[item[0]][1]))
        
    rows = []
    idx_group = 0
    for quantized_int, points in sorted_groups:
        x, y = quantized_to_group[quantized_int]
        points = sorted(points, key=lambda p: p[2])
        points = np.array(points)
        
        # Extract data
        z_coords = points[:, 2].tolist()
        n_ply = len(z_coords)
        
        # Calculate ply thickness: difference between adjacent central points
        # All plies have the same thickness
        if n_ply == 1:
            raise ValueError(f"Single ply laminate is not supported")
        else:
            # Multiple plies: average of differences between adjacent z-coordinates
            z_diffs = [abs(z_coords[i+1] - z_coords[i]) for i in range(n_ply - 1)]
            ply_thickness = np.mean(z_diffs) if len(z_diffs) > 0 else 0.0
            # Total thickness = number of plies * ply thickness
            total_thickness = n_ply * ply_thickness

        rows.append({
            'X': x,
            'Y': y,
            'total_thickness': total_thickness,
            'ply_thickness': ply_thickness,
            'n_ply': n_ply,
            'index': idx_group,
            'Z-list': z_coords,
            'S11-list': points[:, 3].tolist(),
            'S22-list': points[:, 4].tolist(),
            'S33-list': points[:, 5].tolist(),
            'S12-list': points[:, 6].tolist(),
            'S13-list': points[:, 7].tolist(),
            'S23-list': points[:, 8].tolist(),
        })
        idx_group += 1
    
    df = pd.DataFrame(rows)
    
    # Sort by X, then Y
    df = df.sort_values(['X', 'Y']).reset_index(drop=True)
    
    if fname_save is not None:
        df.to_csv(fname_save, index=False)
    
    return df


def check_data_consistency(data_C3D8R: pd.DataFrame, data_SC8R: pd.DataFrame,
            data_S4R: pd.DataFrame, n_samples: int = None) -> bool:
    '''
    Check the data consistency.
    
    Check whether the data points have the same x,y-coordinates,
    and whether the z-coordinate lists are the same.
    
    If n_samples is None, check all data points.
    Otherwise, randomly select n_samples points with the same index from each data source.
    '''
    # Check if all DataFrames have the 'index' column
    for name, df in [('C3D8R', data_C3D8R), ('SC8R', data_SC8R), ('S4R', data_S4R)]:
        if 'index' not in df.columns:
            print(f"Warning: DataFrame {name} does not have 'index' column")
            return False
    
    # Find common indices across all three DataFrames
    indices_C3D8R = set(data_C3D8R['index'].dropna().unique())
    indices_SC8R = set(data_SC8R['index'].dropna().unique())
    indices_S4R = set(data_S4R['index'].dropna().unique())
    
    common_indices = indices_C3D8R & indices_SC8R & indices_S4R
    n_total_indices = len(common_indices)
    
    if n_total_indices == 0:
        print("Error: No common indices found across all three data sources")
        return False
    
    # Select indices (or all if n_samples is None)
    if n_samples is None:
        n_samples = n_total_indices
        selected_indices = list(common_indices)
    else:
        selected_indices = np.random.choice(list(common_indices), n_samples, replace=False)
    
    tolerance = 1e-5
    all_consistent = True
    inconsistent_count = 0
    
    for idx in selected_indices:
        # Get rows with this index from each DataFrame
        row_C3D8R = data_C3D8R[data_C3D8R['index'] == idx]
        row_SC8R = data_SC8R[data_SC8R['index'] == idx]
        row_S4R = data_S4R[data_S4R['index'] == idx]
        
        # Check if exactly one row exists for each index
        if len(row_C3D8R) != 1 or len(row_SC8R) != 1 or len(row_S4R) != 1:
            print(f"Warning: Index {idx} has multiple or no rows in some DataFrames")
            inconsistent_count += 1
            all_consistent = False
            continue
        
        row_C3D8R = row_C3D8R.iloc[0]
        row_SC8R = row_SC8R.iloc[0]
        row_S4R = row_S4R.iloc[0]
        
        # Check X, Y coordinates
        x_C3D8R, y_C3D8R = row_C3D8R['X'], row_C3D8R['Y']
        x_SC8R, y_SC8R = row_SC8R['X'], row_SC8R['Y']
        x_S4R, y_S4R = row_S4R['X'], row_S4R['Y']
        
        coords_match = (np.allclose([x_C3D8R, y_C3D8R], [x_SC8R, y_SC8R], rtol=tolerance) and
                       np.allclose([x_C3D8R, y_C3D8R], [x_S4R, y_S4R], rtol=tolerance))
        
        if not coords_match:
            print(f"Inconsistency at index {idx}: X,Y coordinates don't match")
            print(f"  C3D8R: ({x_C3D8R}, {y_C3D8R})")
            print(f"  SC8R:  ({x_SC8R}, {y_SC8R})")
            print(f"  S4R:   ({x_S4R}, {y_S4R})")
            inconsistent_count += 1
            all_consistent = False
            continue
        
        # Check Z-list
        z_list_C3D8R = row_C3D8R['Z-list']
        z_list_SC8R = row_SC8R['Z-list']
        z_list_S4R = row_S4R['Z-list']
        
        # Convert to numpy arrays for comparison
        z_C3D8R = np.array(z_list_C3D8R)
        z_SC8R = np.array(z_list_SC8R)
        z_S4R = np.array(z_list_S4R)
        
        z_match = (len(z_C3D8R) == len(z_SC8R) == len(z_S4R) and
                  np.allclose(z_C3D8R, z_SC8R, rtol=tolerance) and
                  np.allclose(z_C3D8R, z_S4R, rtol=tolerance))
        
        if not z_match:
            print(f"Inconsistency at index {idx}: Z-list doesn't match")
            print(f"  C3D8R: {z_list_C3D8R}")
            print(f"  SC8R:  {z_list_SC8R}")
            print(f"  S4R:   {z_list_S4R}")
            inconsistent_count += 1
            all_consistent = False
    
    # Print summary
    print(f"\nData consistency check summary:")
    print(f"  Total common indices: {n_total_indices}")
    print(f"  Samples checked: {n_samples}")
    print(f"  Inconsistent samples: {inconsistent_count}")
    print(f"  Consistent samples: {n_samples - inconsistent_count}")
    
    if all_consistent:
        print("  Result: All checked samples are consistent")
    else:
        print("  Result: Some inconsistencies found")
    
    return all_consistent


def save_index_map(data: pd.DataFrame, fname_save: str = 'index-map.dat') -> None:
    '''
    Save the index map of the data points in Tecplot format,
    i.e., 'Variables= X Y index'
    
    Parameters
    ----------
    data: pd.DataFrame
        the data points: 'X', 'Y', 'Z', 'S11', 'S22', 'S33', 'S12', 'S13', 'S23'
    fname_save: str
        the file name to save the index map
    '''
    with open(fname_save, 'w') as f:
        f.write('Variables= X Y index\n')
        for i, row in enumerate(data.values):
            f.write(f'{row[0]:.6f} {row[1]:.6f} {i}\n')


def main():
    '''
    Main function.
    '''
    fname_save_C3D8R = os.path.join(path, 'data', 'stress-field-C3D8R.csv')
    fname_save_SC8R = os.path.join(path, 'data', 'stress-field-SC8R.csv')
    fname_save_S4R = os.path.join(path, 'data', 'stress-field-S4R.csv')
    fname_save_index_map = os.path.join(path, 'data', 'index-map.dat')
    
    print(f"Reading data...")
    data_C3D8R = read_tecplot_data(fname_C3D8R)
    print(f"  C3D8R: {data_C3D8R.shape}")
    data_SC8R = read_tecplot_data(fname_SC8R)
    print(f"  SC8R: {data_SC8R.shape}")
    data_S4R = read_tecplot_data(fname_S4R)
    print(f"  S4R: {data_S4R.shape}")
    
    print(f"Sorting data...")
    data_C3D8R = sort_data_points(data_C3D8R)
    data_SC8R = sort_data_points(data_SC8R)
    data_S4R = sort_data_points(data_S4R)
    
    print(f"Saving data...")
    data_C3D8R = save_as_data_frame(data_C3D8R, fname_save=fname_save_C3D8R)
    data_SC8R = save_as_data_frame(data_SC8R, fname_save=fname_save_SC8R)
    data_S4R = save_as_data_frame(data_S4R, fname_save=fname_save_S4R)
    
    print(f"Checking data consistency...")
    check_data_consistency(data_C3D8R, data_SC8R, data_S4R, n_samples=None)
    
    print(f"Saving index map...")
    save_index_map(data_C3D8R, fname_save=fname_save_index_map)
    
    print(f"Done.")


if __name__ == '__main__':
    main()
