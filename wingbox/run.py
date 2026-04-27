
import os
import platform
import time
import json

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files


N_CPU = 4
FLAG_UVARM = False


if __name__ == '__main__':
    
    t0 = time.time()
    
    clean_pyc_files()

    fname = 'default-parameters.json'
    with open(fname, 'r') as f:
        parameters = json.load(f)

    if platform.system() == 'Windows':
        os.system('abaqus cae script=wingbox_model.py')
    else:
        os.system('abaqus cae noGUI=wingbox_model.py')
    clean_temporary_files()
    
    if parameters['pMesh']['failure_model'] == 'LaRC05':
        print('>>> Running job with LaRC05 failure model...')
        name_job = str(parameters['name_job'])
        os.system('abaqus interactive job=%s user=uvarm.f90 cpus=%d'%(name_job, N_CPU))
        clean_temporary_files()

    t2 = time.time()
    
    print('>>> =============================================')
    print('>>> Time [total]: %.2f min'%((t2-t0)/60.0))
    print('>>> =============================================')

    

