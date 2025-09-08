
import os
import time

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files


N_CPU = 4


if __name__ == '__main__':
    
    t0 = time.time()
    
    clean_pyc_files()

    print('>>> =============================================')
    t1 = time.time()
    
    name_job = 'Job_OHT'
    
    os.system('abaqus cae noGUI=open-hole-test.py')
    
    clean_temporary_files()
    
    os.system('abaqus interactive job=%s user=uvarm.f90 cpus=%d'%(name_job, N_CPU))
    
    clean_temporary_files()
    
    os.system('abaqus cae noGUI=oht-extract-hole-face.py')
    
    clean_temporary_files()

    t2 = time.time()
    
    print('>>> =============================================')
    print('>>> Time [total]: %.2f min'%((t2-t0)/60.0))
    print('>>> =============================================')

    

