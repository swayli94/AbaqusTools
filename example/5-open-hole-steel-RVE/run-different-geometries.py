'''
Run the PBC-3D analysis for different geometries.
'''

import os
import time
import numpy as np
import json
import copy

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files
from AbaqusTools.pbc import PBC_3DOrthotropic


COMMAND = 'abaqus cae noGUI='

fname_py = 'job-pbc-3d.py'


def run_job(index_run: int, parameters: dict) -> dict:
    
    t1 = time.time()

    StiffMatrix = np.zeros([6,6])

    for i in range(6):
        
        parameters['index_strain_vector'] = i
        parameters['index_run'] = index_run
        with open('parameters.json', 'w') as f:
            json.dump(parameters, f, indent=4)
        
        os.system(COMMAND+fname_py)
        
        clean_temporary_files('%d'%(i))
        
        name_job = 'Job_OHP_%d_%d'%(index_run, i)

        with open(name_job+'-RF.dat', 'r') as f:
            lines = f.readlines()
            
            scale = float(lines[12+i].split()[1])
            for j in range(6):
                StiffMatrix[j,i] = float(lines[j].split()[1])/scale

    results = PBC_3DOrthotropic.calculate_engineering_constants(StiffMatrix)
    
    volume_box = parameters['pGeo']['len_x_plate'] * parameters['pGeo']['len_y_plate'] * parameters['pGeo']['len_z_plate']
    volume_hole = np.pi*parameters['pGeo']['r_hole']**2*parameters['pGeo']['len_z_plate']
    volume = volume_box - volume_hole
    
    results['volume'] = volume
    results['volume_box'] = volume_box
    results['volume_hole'] = volume_hole
    results['run_time (min)'] = (time.time()-t1)/60.0

    print('>>> Job %d: %.2f min'%(index_run+1, results['run_time (min)']))

    return results


if __name__ == '__main__':
    
    t0 = time.time()
    
    clean_pyc_files()
    
    with open('parameters.json', 'r') as f:
        default_parameters = json.load(f)


    index_run = 0
    all_results = {}
    
    for len_x_plate in [50, 100, 200]:
        for len_y_plate in [50, 100]:
            
            area = len_x_plate * len_y_plate
            plate_seedPart_size = np.sqrt(area) / 10.0
            
            for len_z_plate, num_element_thickness in [(5, 3), (10, 5), (50, 10)]:
                
                for r_hole in [5, 10, 20]:
                    
                    parameters = copy.deepcopy(default_parameters)
                    
                    if len_x_plate <= 4*r_hole or len_y_plate <= 4*r_hole:
                        continue
                    
                    parameters['pGeo']['len_x_plate'] = len_x_plate
                    parameters['pGeo']['len_y_plate'] = len_y_plate
                    parameters['pGeo']['len_z_plate'] = len_z_plate
                    parameters['pGeo']['r_hole'] = r_hole
                    
                    parameters['pMesh']['num_element_thickness'] = num_element_thickness
                    parameters['pMesh']['plate_seedPart_size'] = plate_seedPart_size
                    
                    results = run_job(index_run, parameters)
                    results.update(parameters)
                    
                    all_results[index_run] = copy.deepcopy(results)

                    index_run += 1

                    with open('all-results.json', 'w') as f:
                        json.dump(all_results, f, indent=4)

    print('>>> =============================================')
    print('>>> Time [total]: %.2f min'%((time.time()-t0)/60.0))
    print('>>> =============================================')

    



