
import os
import platform
import sys
import time

from AbaqusTools.functions import clean_pyc_files, clean_temporary_files

from params import get_parameter_file, load_parameters, get_failure_model


def submits_job_in_cae(parameters):
    '''
    Whether `wingbox_model.py` submits the analysis job from inside Abaqus/CAE.

    The temporary-file cleanup removes `*.com`, `*.env`, `*.sim` and `*.stt`,
    which the solver still needs, so it must never run while a job is alive.
    '''
    if parameters.get('not_run_job', False):
        return False
    if get_failure_model(parameters['pMesh']) == 'larc05':
        return False
    return str(parameters.get('execution_mode', 'default')).lower() in ('default', 'auto')


def get_postprocess_parameters(parameters):
    '''
    Get the failure post-processing settings, see `postprocess_failure.py`
    (mode 'full') and `extract_results.py` (mode 'fast').
    '''
    postprocess = parameters.get('postprocess', {})
    return {
        'enabled': bool(postprocess.get('enabled', True)),
        'mode': str(postprocess.get('mode', 'full')).lower(),
        'slim_odb': bool(postprocess.get('slim_odb', True)),
        'delete_source_odb': bool(postprocess.get('delete_source_odb', False)),
    }


if __name__ == '__main__':

    t0 = time.time()

    clean_pyc_files()

    fname = get_parameter_file(sys.argv)
    print('>>> Reading parameters from: %s' % fname)
    parameters = load_parameters(fname)

    name_job = str(parameters['name_job'])

    #* A stale output database of a previous run would otherwise be reduced and
    #* reported as if it belonged to this design.
    for suffix in ('.odb', '_failure_envelope.npz', '_failure_summary.json'):
        if os.path.isfile(name_job + suffix):
            os.remove(name_job + suffix)

    if platform.system() == 'Windows':
        command = 'abaqus cae script=wingbox_model.py -- --params %s' % fname
    else:
        command = 'abaqus cae noGUI=wingbox_model.py -- --params %s' % fname
    status = os.system(command)
    if status != 0:
        raise RuntimeError('Abaqus/CAE failed (exit status %d): %s' % (status, command))

    #* `wingbox_model.py` waits for the job by default, so the analysis is over
    #* once Abaqus/CAE returns.  When waiting is switched off the job may still
    #* be running, and cleaning up now would delete the files it is using.
    if submits_job_in_cae(parameters) and not parameters.get('wait_for_completion', True):
        print('>>> --------------------')
        print('    [Skip] clean_temporary_files()')
        print('    The job was submitted with "wait_for_completion": false,')
        print('    so *.com/*.env/*.sim/*.stt are left in place for the solver.')
        print('>>>')
    else:
        clean_temporary_files()

    if get_failure_model(parameters['pMesh']) == 'larc05':
        # The LaRC05 deck is patched after the model is built, so the job is
        # submitted here instead of from inside Abaqus/CAE.
        # standard_parallel=solver keeps the element loop serial: the LaRC05
        # user subroutine is not thread-safe (module-level state), while the
        # linear solver still uses all cpus.
        print('>>> Running job with LaRC05 failure model...')
        #* The solver runs in the scratch directory; tell module larc05Track
        #* in uvarm.f90 where to write larc05_fi_track_<pid>.txt.
        os.environ['LARC05_TRACK_DIR'] = os.path.abspath('.')
        command = ('abaqus interactive job=%s user=uvarm.f90 cpus=%d '
                   'standard_parallel=solver'
                   % (name_job, parameters['pRun']['numCpus']))
        status = os.system(command)
        clean_temporary_files()
        if status != 0:
            raise RuntimeError('Abaqus job failed (exit status %d): %s' % (status, command))

    t1 = time.time()

    #* Reduce the analysis output to the design-evaluation numbers
    postprocess = get_postprocess_parameters(parameters)
    if postprocess['enabled'] and os.path.isfile(name_job + '.odb'):
        if postprocess['mode'] == 'fast':
            #* max FI from the UVARM tracking file, eigenvalues from the
            #* data file, tip displacement from the nodal U field.
            command = 'abaqus python extract_results.py --job %s' % name_job
        else:
            command = 'abaqus python postprocess_failure.py --job %s' % name_job
            if postprocess['slim_odb']:
                command += ' --slim-odb'
            if postprocess['delete_source_odb']:
                command += ' --delete-source-odb'
        print('>>> %s' % command)
        status = os.system(command)
        if status != 0:
            raise RuntimeError(
                'Failure post-processing failed (exit status %d): %s' % (status, command))

    t2 = time.time()

    print('>>> =============================================')
    print('>>> Time [analysis]:      %.2f min' % ((t1-t0)/60.0))
    print('>>> Time [postprocess]:   %.2f min' % ((t2-t1)/60.0))
    print('>>> Time [total]:         %.2f min' % ((t2-t0)/60.0))
    print('>>> =============================================')
