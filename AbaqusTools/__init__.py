'''
You are suggested to install `abqpy` to provide type hints for Abaqus/Python scripting.

    `abqpy` is a Python package providing type hints for Python scripting of Abaqus, 
    you can use it to write your Python script of Abaqus fluently, even without doing anything in Abaqus. 
    It also provides some simple APIs to execute the Abaqus commands so that you can run your 
    Python script to build the model, submit the job and extract the output data in just one Python script, 
    even without opening the Abaqus/CAE.

    https://github.com/haiiliin/abqpy
'''

#* Only import Abaqus modules when it is called by Abaqus
import sys

cmd_arguments = str(sys.argv)

if '-tmpdir' in cmd_arguments or '-cae' in cmd_arguments or '-interactive' in cmd_arguments:

    IS_ABAQUS = True
    
else:
    
    IS_ABAQUS = False
    
if IS_ABAQUS:
    
    from AbaqusTools.part import Part
    from AbaqusTools.model import Model, NodeOperation
    from AbaqusTools.odb import OdbOperation
    from AbaqusTools.functions import load_parameters, clean_pyc_files, clean_temporary_files
    