'''
Run the analysis for Open Hole Plate (OHP) specimen.
'''

import os
import time
import json

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files
from implicit_modelling import update_parameters

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
        
        os.system(COMMAND+fname_py)
        
        clean_temporary_files('%d'%(i))
        
        name_job = 'Job_OHP_%d_%d'%(index_run, i)

        t2 = time.time()
        
        print('>>> Time [strain vector %d]: %.2f min'%(i, (t2-t1)/60.0))
        print('>>> =============================================')

