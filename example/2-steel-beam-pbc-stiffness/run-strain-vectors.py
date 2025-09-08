'''
Run four load cases with the basis strain vectors

>>> cp -f ./example/2-beam-pbc-steel/*.py ./ & clean.bat & python run-strain-vectors.py
'''

import os
import time
import numpy as np

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files


COMMAND = 'abaqus cae noGUI='

fname_py = 'beam-strain-vector-C3D8R.py'


if __name__ == '__main__':
    
    t0 = time.time()
    
    clean_pyc_files()

    StiffMatrix = np.zeros([4,4])

    for i in range(4):
        
        print('>>> =============================================')
        t1 = time.time()
        
        with open('temp-strain-vector.txt', 'w') as f:
            f.write('  %d \n'%(i))
        
        os.system(COMMAND+fname_py)
        
        clean_temporary_files('%d'%(i))
        
        name_job = 'Job_beam_%d'%(i)

        with open(name_job+'-RF.dat', 'r') as f:
            lines = f.readlines()
            
            scale = float(lines[10+i].split()[1])
            for j in range(4):
                StiffMatrix[j,i] = float(lines[j].split()[1])/scale
        
        t2 = time.time()
        
        print('>>> Time [strain vector %d]: %.2f min'%(i, (t2-t1)/60.0))
        print('>>> =============================================')


    with open('stiffness-matrix.dat', 'w') as f:
        for j in range(4):
            for i in range(4):
                f.write('%15.3E'%(StiffMatrix[j,i]))
            f.write('\n')

    print('>>> =============================================')
    print('>>> Time [total]: %.2f min'%((t2-t0)/60.0))
    print('>>> =============================================')

    

