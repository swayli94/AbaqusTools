import os
import time
import json

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files
from implicit_modelling import update_parameters

COMMAND = 'abaqus cae noGUI='

fname_py = 'mif_single_lap_S4R.py'

DISPLACEMENT = [[0.01, 0.0, 0.0],
                [0.0, 0.01, 0.0],
                [0.0, 0.0, 0.01]]


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

        print('>>> Time [Case %d]: %.2f min'%(i_case, (time.time()-t1)/60.0))

