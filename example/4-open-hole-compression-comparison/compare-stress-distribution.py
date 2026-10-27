'''
Compare the thickness-distributed stress distribution of the three data sources.

Randomly or manually select a few data point locations,
find the corresponding data points in the three data sources,
and compare the stress distribution.

Also, reconstruct the stress distribution from the data of central point of each ply.
For better visualization, we need to plot the stress of each ply with three points,
which have the same x,y-coordinates, different z-coordinates,
and the same stress values as the central point of each ply.
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
import matplotlib.pyplot as plt
import ast

N_COL = 9

fname_C3D8R = 'stress-field-C3D8R.csv'
fname_SC8R = 'stress-field-SC8R.csv'
fname_S4R = 'stress-field-S4R.csv'

fname_C3D8R = os.path.join(path, 'data', fname_C3D8R)
fname_SC8R = os.path.join(path, 'data', fname_SC8R)
fname_S4R = os.path.join(path, 'data', fname_S4R)


def read_data(fname: str) -> pd.DataFrame:
    '''
    Read the data from the CSV file, restore the data as int, float, list, etc.
    '''
    data_frame = pd.read_csv(fname)
    
    # Parse string representations of lists back to actual lists
    # CSV stores lists as strings like "[0.05, 0.15, ...]", need to convert them back
    list_columns = [col for col in data_frame.columns if col.endswith('-list')]
    for col in list_columns:
        data_frame[col] = data_frame[col].apply(ast.literal_eval)
    
    return data_frame


def plot_stress_distribution(ax: plt.Axes, data_frame: pd.DataFrame, idx: int,
                    variable: str='S11',
                    label: str=None, color: str='black',
                    linestyle: str='-', linewidth: float=1.5) -> None:
    '''
    Plot the stress distribution of the data points.
    
    Parameters
    ----------
    ax: plt.Axes
        the axes to plot the stress distribution
    data_frame: pd.DataFrame
        the data frame containing the data points
    idx: int
        the index of the data point
    variable: str
        the variable to plot
    label: str
        the label of the stress distribution
    color: str
        the color of the stress distribution
    linestyle: str
        the linestyle of the stress distribution
    linewidth: float
        the linewidth of the stress distribution
    '''
    Z_list = data_frame['Z-list'].iloc[idx]
    variable_list = data_frame[f'{variable}-list'].iloc[idx]
    ply_thickness = data_frame['ply_thickness'].iloc[idx]
    n_ply = len(Z_list)
    half_ply_thickness = ply_thickness / 2.0
    
    # Construct the stress distribution
    stress_distribution = np.zeros((n_ply*2, 2))
    for i_ply in range(n_ply):
        z_center = Z_list[i_ply]
        stress_distribution[2*i_ply,   1] = z_center - half_ply_thickness
        stress_distribution[2*i_ply+1, 1] = z_center + half_ply_thickness
        stress_distribution[2*i_ply,   0] = variable_list[i_ply]
        stress_distribution[2*i_ply+1, 0] = variable_list[i_ply]
    
    # Plot the stress distribution
    ax.plot(stress_distribution[:, 0], stress_distribution[:, 1],
            linestyle=linestyle, linewidth=linewidth, color=color, label=label)


def plot_stress_comparison(data_C3D8R: pd.DataFrame, data_SC8R: pd.DataFrame,
                        data_S4R: pd.DataFrame, idx: int,
                        fig_title: str='Stress Distribution Comparison',
                        fname_save: str=None) -> None:
    '''
    Plot the stress comparison of the three data sources.
    '''
    fig = plt.figure(figsize=(16, 10))
    
    variable_list = ['S11', 'S22', 'S12', 'S33', 'S13', 'S23']
    
    handles = None
    labels = None
    
    for i_var, variable in enumerate(variable_list):
        ax = fig.add_subplot(2, 3, i_var + 1)
        plot_stress_distribution(ax, data_C3D8R, idx, variable=variable, label='C3D8R', color='black', linewidth=2.0)
        plot_stress_distribution(ax, data_SC8R, idx, variable=variable, label='SC8R', color='orange', linewidth=1.0)
        plot_stress_distribution(ax, data_S4R, idx, variable=variable, label='S4R', color='green', linewidth=1.0, linestyle='--')
        ax.set_xlabel(variable)
        ax.set_ylabel('Z')
        
        if i_var == 0:
            handles, labels = ax.get_legend_handles_labels()
    
    x = data_C3D8R['X'].iloc[idx]
    y = data_C3D8R['Y'].iloc[idx]
    note = f' (index={idx}, x={x:.2f}, y={y:.2f})'
    
    fig.suptitle(fig_title+note, fontsize=14, y=0.95)
    
    fig.legend(handles, labels, 
               loc='upper center', 
               ncol=3, 
               bbox_to_anchor=(0.5, 0.93),
               frameon=True)
    
    if fname_save is not None:
        plt.savefig(fname_save, dpi=100, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def main():
    '''
    Main function.
    '''
    data_C3D8R = read_data(fname_C3D8R)
    data_SC8R = read_data(fname_SC8R)
    data_S4R = read_data(fname_S4R)
    
    os.makedirs(os.path.join(path, 'figure'), exist_ok=True)
    
    index_list = [0, 44, 157, 480, 490, 510, 481, 491, 511, 1023, 979, 866]
    for idx in index_list:
        fname_save = os.path.join(path, 'figure', f'stress-comparison-{idx}.png')
        plot_stress_comparison(data_C3D8R, data_SC8R, data_S4R, idx, fname_save=fname_save)
        

if __name__ == '__main__':
    main()
