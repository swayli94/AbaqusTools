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
SCALE = 1E-6


if __name__ == '__main__':
    
    t0 = time.time()
    
    clean_pyc_files()

    StiffMatrix = np.zeros([6,6])

    for i in range(6):
        
        print('>>> =============================================')
        t1 = time.time()
        
        with open('temp-strain-vector.txt', 'w') as f:
            f.write('  %d \n'%(i))
        
        os.system(COMMAND+fname_py)
        
        clean_temporary_files('%d'%(i))
        
        name_job = 'Job_beam_%d'%(i)

        with open(name_job+'-RF.dat', 'r') as f:
            lines = f.readlines()
            
            for j in range(6):
                StiffMatrix[j,i] = float(lines[j].split()[1])/SCALE
        
        t2 = time.time()
        
        print('>>> Time [strain vector %d]: %.2f min'%(i, (t2-t1)/60.0))
        print('>>> =============================================')


    engineering_constants = PBC_3DOrthotropic.calculate_engineering_constants(StiffMatrix)

    with open('homogenized-properties.json', 'w') as f:
        json.dump(engineering_constants, f, indent=4)

    print('>>> =============================================')
    print('>>> Time [total]: %.2f min'%((t2-t0)/60.0))
    print('>>> =============================================')

    

