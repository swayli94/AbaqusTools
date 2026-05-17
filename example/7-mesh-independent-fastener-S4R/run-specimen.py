import os
import time
import json
import numpy as np

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files
from implicit_modelling import (update_parameters, read_tecplot,
                read_ctf_from_rf_dat, calculate_bypass_load,
                calculate_ply_level_stress_field, save_tecplot)

COMMAND = 'abaqus cae noGUI='

fname_py = 'mif_single_lap_S4R.py'

DISPLACEMENT = [[0.01, 0.0, 0.0],
                [0.0, 0.01, 0.0],
                [0.0, 0.0, 0.01]]

NAME_INSTANCES = ['PLATE_0', 'PLATE_1']
NAME_SET = 'PARTITION_SQUARE'


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

    summary = {}

    for i_case in range(3):
        
        t1 = time.time()
        
        default_parameters['index_run'] = 0
        default_parameters['index_case'] = i_case
        default_parameters['displacement'] = DISPLACEMENT[i_case]
        update_parameters(default_parameters, target_volume_fraction,
                        use_implicit_modelling=False)
        with open('parameters.json', 'w') as f:
            json.dump(default_parameters, f, indent=4)
                    
        os.system(COMMAND+fname_py)
        
        clean_temporary_files('%d'%(i_case))

        # ---------------------
        # Post-process
        # ---------------------
        pMesh = default_parameters['pMesh']
        pFastener = default_parameters['pGeo']['fasteners'][0]
        r_hole = pFastener['r_hole']
        thickness = default_parameters['pGeo']['len_z_plate']
        characteristic_distance = r_hole * 2.0

        name_job = 'Job_MIF_%d_%d'%(index_run, i_case)
        fname_midplane = name_job + '-mid-plane.dat'
        fname_rf = name_job + '-RF.dat'
        if not os.path.exists(fname_midplane) or not os.path.exists(fname_rf):
            print(f'File not found: {fname_midplane} or {fname_rf}')
            print('>>> Time [Case %d]: %.2f min'%(i_case, (time.time()-t1)/60.0))
            continue
        
        #* Read bearing load from RF data
        ctf_x, ctf_y, ctf_z = read_ctf_from_rf_dat(fname_rf, fastener_index=0)
        bearing_force = (ctf_x**2 + ctf_y**2)**0.5
        angle_degree = np.rad2deg(np.arctan2(ctf_y, ctf_x))
        
        summary[i_case] = {
            'bearing_force': bearing_force,
            'angle_degree': angle_degree,
            'ctf_x': ctf_x,
            'ctf_y': ctf_y,
            'ctf_z': ctf_z,
        }
        
        #* Load mid-plane data for partition square
        zones = read_tecplot(fname_midplane)
        for i_zone, zone in enumerate(zones):
   
            #* Calculate bypass load
            results = calculate_bypass_load(
                X=zone['X'], Y=zone['Y'],
                N11=zone['N11'], N22=zone['N22'], N12=zone['N12'],
                bearing_force=bearing_force,
                angle_degree=angle_degree
            )
            
            summary[i_case]['zone_%d'%i_zone] = {
                'N_x_bypass': results['N_x_bypass'],
                'N_y_bypass': results['N_y_bypass'],
                'N_xy_bypass': results['N_xy_bypass'],
            }
            
            #* Calculate ply-level stress field
            field, results_by_plies = calculate_ply_level_stress_field(
                pMesh=pMesh,
                r_hole=r_hole,
                total_thickness=thickness,
                N11=results['N_x_bypass'],
                N22=results['N_y_bypass'],
                N12=results['N_xy_bypass'],
                bearing_load=bearing_force,
                angle_load_degree=angle_degree,
                characteristic_distance=characteristic_distance,
                n_points_radial=32,
                n_points_angular=64,
            )

            save_tecplot(name_job+'-stress-field-%d.dat'%(i_zone),
                    field, results_by_plies)

        run_time = time.time() - t1
        summary[i_case]['time'] = run_time

        print('>>> Time [Case %d]: %.2f min'%(i_case, run_time/60.0))

    with open('summary-%d.json'%(index_run), 'w') as f:
        json.dump(summary, f, indent=4)
    