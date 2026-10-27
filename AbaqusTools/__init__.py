'''
Type hints for Abaqus/Python scripting come from `abqpy`, vendored as a git
submodule in `external/abqpy` and pinned to its 2023 branch to match Abaqus 2023.

`pyrightconfig.json` puts `external/abqpy/src` on the type checker's path, which
is all that is needed for completion and signatures in the editor.

Do NOT `pip install abqpy` for this package and do NOT import its modules at
runtime: `abqpy`'s `abaqus` module calls `abqpy.run(cae=True)` when imported,
which re-launches the running script through `abaqus cae`. The `IS_ABAQUS` guard
below is what keeps that from happening in ordinary Python.

    https://github.com/haiiliin/abqpy
'''

#* Only import Abaqus modules when it is called by Abaqus
import sys

cmd_arguments = str(sys.argv)

if '-tmpdir' in cmd_arguments or '-cae' in cmd_arguments or '-interactive' in cmd_arguments:

    IS_ABAQUS = True
    
else:
    
    IS_ABAQUS = False
    
#* Pure data, no Abaqus import, usable in a normal Python as well
from AbaqusTools.materials import MATERIAL_LIBRARY, get_material, list_materials

if IS_ABAQUS:

    from AbaqusTools.part import Part
    from AbaqusTools.model import Model, NodeOperation
    from AbaqusTools.odb import OdbOperation
    from AbaqusTools.functions import load_parameters, clean_pyc_files, clean_temporary_files
    