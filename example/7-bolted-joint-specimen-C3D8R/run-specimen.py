'''
Run the analysis for Open Hole Plate (OHP) specimen.
'''

import os
import time
import json

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files

COMMAND = 'abaqus cae noGUI='

fname_py = 'bolted_joint_single_lap_C3D8R.py'

DISPLACEMENT = [[1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0]]


if __name__ == '__main__':
    
    t0 = time.time()
    
    clean_pyc_files()

    with open('default-parameters.json', 'r') as f:
        default_parameters = json.load(f)
    
    index_run = default_parameters['index_run']

    for i_case in range(3):
        
        t1 = time.time()
        
        default_parameters['index_case'] = i_case
        default_parameters['displacement'] = DISPLACEMENT[i_case]
        with open('parameters.json', 'w') as f:
            json.dump(default_parameters, f, indent=4)
            
        os.system(COMMAND+fname_py)
        
        clean_temporary_files('%d'%(i_case))
        
        name_job = 'Job_OHP_%d_%d'%(index_run, i_case)

        t2 = time.time()
        
        print('>>> Time [strain vector %d]: %.2f min'%(i_case, (t2-t1)/60.0))
