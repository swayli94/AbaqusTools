'''
Parameter-file helpers shared by the driver script and the model script.

`run.py` runs in the system Python, `wingbox_model.py` runs in the Abaqus
Python interpreter, and both need to agree on which JSON file is being used.
Keep this module free of numpy/scipy and Abaqus imports so that either
interpreter can import it.
'''
import json


DEFAULT_PARAMETER_FILE = 'parameters-simple-wingbox.json'

#* Failure analysis is optional: `pMesh['failure_model']` can be absent, null,
#* or 'none', in which case no failure criterion is defined and no failure
#* index is written to the output database.
FAILURE_MODELS = ('hashin', 'larc05', 'none')


def get_failure_model(pMesh):
    '''
    Normalised failure model name: 'hashin', 'larc05' or 'none'.
    '''
    failure_model = pMesh.get('failure_model')

    if failure_model is None:
        return 'none'

    failure_model = str(failure_model).strip().lower()
    if failure_model not in FAILURE_MODELS:
        raise ValueError('Invalid failure_model: %s, should be one of %s.'
                         % (failure_model, FAILURE_MODELS))

    return failure_model


def get_parameter_file(argv):
    '''
    Get parameter JSON path from command-line arguments.

    Default:
        parameters-simple-wingbox.json

    Optional:
        --params path/to/parameters.json

    Arguments after a bare `--` are the ones forwarded by
    `abaqus cae noGUI=script.py -- --params file.json`.
    '''
    if '--' in argv:
        argv = argv[argv.index('--') + 1:]

    if '--params' in argv:
        index = argv.index('--params')
        if index + 1 >= len(argv):
            raise ValueError('Missing value after --params.')
        return argv[index + 1]

    return DEFAULT_PARAMETER_FILE


def load_parameters(fname):
    '''
    Read a parameter JSON file.
    '''
    with open(fname, 'r') as f:
        return json.load(f)
