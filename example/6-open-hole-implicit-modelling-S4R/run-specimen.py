'''
Run the analysis for Open Hole Plate (OHP) specimen.
'''

import os
import time
import json

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files
from implicit_modelling import (update_parameters,
                    calculate_ply_level_stress_field, save_tecplot)

COMMAND = 'abaqus cae noGUI='

fname_py = 'job-specimen.py'


if __name__ == '__main__':
    
    t0 = time.time()
    
    clean_pyc_files()

    with open('default-parameters.json', 'r') as f:
        default_parameters = json.load(f)
    
    index_run = default_parameters['index_run']
    if 'target_volume_fraction' in default_parameters:
        target_volume_fraction=default_parameters['target_volume_fraction']
    else:
        target_volume_fraction=0.4

    for i in range(3):
        
        print('>>> =============================================')
        t1 = time.time()
        
        default_parameters['index_strain_vector'] = i
        update_parameters(default_parameters, target_volume_fraction)
        with open('parameters.json', 'w') as f:
            json.dump(default_parameters, f, indent=4)
            
        scale = default_parameters['strain_scale']
        # characteristic_distance = default_parameters['characteristic_distance']
        characteristic_distance = default_parameters['pGeo']['r_hole']*2.0
        
        os.system(COMMAND+fname_py)
        
        clean_temporary_files('%d'%(i))
        
        # -----------------------------------
        # Read section stress
        # -----------------------------------
        name_job = 'Job_OHP_%d_%d'%(index_run, i)
        fname = name_job + '-mid-plane.dat'
        if not os.path.exists(fname):
            print(f'File not found: {fname}')
            print('>>> Time [Case %d]: %.2f min'%(i, (time.time()-t1)/60.0))
            continue
            
        with open(fname, 'r') as f:
            lines = f.readlines()[2:]
    
            for k, line in enumerate(lines):
                values = line.strip().split()
                if len(values) == 0:
                    continue

                field, results_by_plies = calculate_ply_level_stress_field(
                    parameters=default_parameters,
                    N11=float(values[4]),
                    N22=float(values[5]),
                    N12=float(values[6]),
                    characteristic_distance=characteristic_distance,
                    n_points_radial=32,
                    n_points_angular=64,
                )

                save_tecplot(name_job+'-stress-field.dat', field, results_by_plies)

        print('>>> Time [Case %d]: %.2f min'%(i, (time.time()-t1)/60.0))

