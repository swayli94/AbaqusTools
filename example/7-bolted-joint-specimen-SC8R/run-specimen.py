'''
Run the analysis for Open Hole Plate (BJ) specimen.
'''

import os
import time
import json

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files

COMMAND = 'abaqus cae noGUI='

fname_py = 'bolted_joint_single_lap_SC8R.py'

DISPLACEMENT = [[0.01, 0.0, 0.0],
                [0.0, 0.01, 0.0],
                [0.0, 0.0, 0.01]]


if __name__ == '__main__':
    
    t0 = time.time()
    
    clean_pyc_files()

    with open('default-parameters.json', 'r') as f:
        default_parameters = json.load(f)
    
    for i_case in range(3):
        
        t1 = time.time()
        
        default_parameters['index_run'] = 0
        default_parameters['index_case'] = i_case
        default_parameters['displacement'] = DISPLACEMENT[i_case]
        with open('parameters.json', 'w') as f:
            json.dump(default_parameters, f, indent=4)
            
        os.system(COMMAND+fname_py)
        
        clean_temporary_files('%d'%(i_case))
        
        t2 = time.time()
        
        print('>>> Time [strain vector %d]: %.2f min'%(i_case, (t2-t1)/60.0))
