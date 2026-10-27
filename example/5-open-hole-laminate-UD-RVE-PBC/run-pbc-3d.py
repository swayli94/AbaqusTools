'''
Run four load cases with the basis strain vectors

>>> cp -f ./example/2-beam-pbc-steel/*.py ./ & clean.bat & python run-strain-vectors.py
'''

import os
import time
import numpy as np
import json

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files
from AbaqusTools.pbc import PBC_3DOrthotropic


COMMAND = 'abaqus cae noGUI='

fname_py = 'job-pbc-3d.py'


if __name__ == '__main__':
    
    t0 = time.time()
    
    clean_pyc_files()

    with open('default-parameters.json', 'r') as f:
        default_parameters = json.load(f)
    
    index_run = default_parameters['index_run']
    StiffMatrix = np.zeros([6,6])

    for i in range(6):
        
        print('>>> =============================================')
        t1 = time.time()
        
        default_parameters['index_strain_vector'] = i
        with open('parameters.json', 'w') as f:
            json.dump(default_parameters, f, indent=4)
            
        scale = default_parameters['strain_scale']
        
        os.system(COMMAND+fname_py)
        
        clean_temporary_files('%d'%(i))
        
        name_job = 'Job_OHP_%d_%d'%(index_run, i)

        with open(name_job+'-RF.dat', 'r') as f:
            lines = f.readlines()
            
            for j in range(6):
                StiffMatrix[j,i] = float(lines[j].split()[1])/scale
        
        t2 = time.time()
        
        print('>>> Time [strain vector %d]: %.2f min'%(i, (t2-t1)/60.0))
        print('>>> =============================================')


    engineering_constants = PBC_3DOrthotropic.calculate_engineering_constants(StiffMatrix)

    with open('homogenized-properties.json', 'w') as f:
        json.dump(engineering_constants, f, indent=4)

    print('>>> =============================================')
    print('>>> Time [total]: %.2f min'%((t2-t0)/60.0))
    print('>>> =============================================')

    

