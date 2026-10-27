
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
    
    os.system('abaqus cae noGUI=open-hole-compression.py')
    
    clean_temporary_files()
    
    # os.system('abaqus interactive job=%s user=uvarm.f90 cpus=%d'%(name_job, N_CPU)) # failure_model = LaRC05, user_subroutine = UVARM
    os.system('abaqus interactive job=%s cpus=%d'%(name_job, N_CPU))
    
    clean_temporary_files()
    
    os.system('abaqus cae noGUI=extract-data-C3D8R.py')
    
    clean_temporary_files()

    t2 = time.time()
    
    print('>>> =============================================')
    print('>>> Time [total]: %.2f min'%((t2-t0)/60.0))
    print('>>> =============================================')

    

