'''
Reduce a wingbox analysis output database to what a failure analysis needs.

The failure indices are evaluated at every section point of every ply, so the
raw output database grows with the number of plies (a 100-ply cover writes 300
values per element with a three-point Simpson rule).  For failure analysis only
the through-thickness envelope matters, so this script

- reads the failure indices of the last static frame,
- reduces them to one value per element, `max(|value|)` over all section points,
- does the same for the von Mises stress when a stress field is present, i.e.,
  for the metallic ribs, which carry no failure index,
- writes them to a compact `.npz` file plus a `.json` summary, and
- optionally rebuilds a slim output database that keeps only the first buckling
  mode shape, the deformed shape of the static step, and the envelope fields.

The slim database is typically two orders of magnitude smaller than the source,
which can then be deleted.

Usage (Abaqus Python, no CAE license needed)
--------------------------------------------
    abaqus python postprocess_failure.py --job Job_WB [--slim-odb]
                                         [--delete-source-odb]
'''
import os
import sys
import json

import numpy as np

import odbAccess
from odbAccess import openOdb
from abaqusConstants import (THREE_D, DEFORMABLE_BODY, TIME, NODAL, SCALAR,
                             VECTOR, MAGNITUDE, MISES, WHOLE_ELEMENT, CENTROID)


#* Failure indices of the built-in Hashin model and of the LaRC05 subroutine.
#* `UVARM` fields are named UVARM1, UVARM2, ... in the output database.
HASHIN_VARIABLES = ('DMICRT', 'HSNFTCRT', 'HSNFCCRT', 'HSNMTCRT', 'HSNMCCRT')


def parse_arguments(argv):
    '''
    Minimal command-line parsing (`argparse` is awkward in Abaqus Python 2.7).
    '''
    arguments = {
        'job': None,
        'odb': None,
        'slim_odb': False,
        'delete_source_odb': False,
    }
    i = 1
    while i < len(argv):
        item = argv[i]
        if item == '--job':
            arguments['job'] = argv[i+1]; i += 2
        elif item == '--odb':
            arguments['odb'] = argv[i+1]; i += 2
        elif item == '--slim-odb':
            arguments['slim_odb'] = True; i += 1
        elif item == '--delete-source-odb':
            arguments['delete_source_odb'] = True; i += 1
        elif item == '--':
            i += 1
        else:
            raise ValueError('Unknown argument: %s' % item)

    if arguments['odb'] is None:
        if arguments['job'] is None:
            raise ValueError('Either --job or --odb is required.')
        arguments['odb'] = arguments['job'] + '.odb'

    if arguments['job'] is None:
        arguments['job'] = os.path.splitext(os.path.basename(arguments['odb']))[0]

    return arguments


def get_static_step(odb):
    '''
    Get the name of the last general static step, i.e., the loaded state.
    '''
    name_static = None
    for name in odb.steps.keys():
        if 'BUCKLE' not in str(odb.steps[name].procedure).upper():
            name_static = name
    return name_static


def get_buckle_step(odb):
    '''
    Get the name of the linear perturbation buckling step.
    '''
    for name in odb.steps.keys():
        if 'BUCKLE' in str(odb.steps[name].procedure).upper():
            return name
    return None


def get_failure_variable_names(frame):
    '''
    Get the failure index fields present in a frame.

    Both the built-in Hashin criteria and the LaRC05 user-defined output
    variables (UVARM1, UVARM2, ...) are recognised.
    '''
    names = []
    for name in frame.fieldOutputs.keys():
        name = str(name)
        if name in HASHIN_VARIABLES or name.startswith('UVARM'):
            names.append(name)
    return sorted(names)


def get_stress_field_name(frame):
    '''
    Name of the von Mises envelope field, when a stress field is available.

    Metallic ribs carry no composite layup, so their strength is checked from
    the stress field instead of a failure index.
    '''
    if 'S' in frame.fieldOutputs.keys():
        return 'S'
    return None


def compute_element_envelope(field):
    '''
    Reduce a section-point field to one value per element.

    Parameters
    ----------------
    field: FieldOutput
        field output of one failure index, evaluated at the section points of
        the composite layups.

    Returns
    ----------------
    envelope: dict
        {instance name: (element labels [n], max|value| [n])}

    Notes
    ----------------
    `bulkDataBlocks` returns the values of one instance / element type /
    section point as numpy arrays, which keeps the reduction vectorised.
    Iterating over `field.values` instead would take hours for a model with
    millions of section-point values.
    '''
    labels_of_instance = {}
    values_of_instance = {}

    for block in field.bulkDataBlocks:

        if block.elementLabels is None or len(block.elementLabels) == 0:
            continue

        name_instance = 'ASSEMBLY'
        if block.instance is not None:
            name_instance = str(block.instance.name)

        labels = np.asarray(block.elementLabels, dtype=np.int64)
        data = np.asarray(block.data, dtype=np.float64)
        if data.ndim == 1:
            data = data.reshape(-1, 1)
        # Several integration points of the same element appear as consecutive
        # entries with a repeated label, which the reduction below takes care of.
        values = np.abs(data).max(axis=1)

        labels_of_instance.setdefault(name_instance, []).append(labels)
        values_of_instance.setdefault(name_instance, []).append(values)

    envelope = {}
    for name_instance in labels_of_instance:
        labels = np.concatenate(labels_of_instance[name_instance])
        values = np.concatenate(values_of_instance[name_instance])
        unique_labels = np.unique(labels)
        reduced = np.zeros(unique_labels.shape, dtype=np.float64)
        np.maximum.at(reduced, np.searchsorted(unique_labels, labels), values)
        envelope[name_instance] = (unique_labels, reduced)

    return envelope


def read_eigenvalues(fname_dat):
    '''
    Read the buckling eigenvalues from the `*.dat` file.

    Only the requested mode shapes are written to the output database, while
    all eigenvalues are always printed to the data file.
    '''
    eigenvalues = []
    if not os.path.isfile(fname_dat):
        return eigenvalues

    in_table = False
    with open(fname_dat, 'r') as f:
        for line in f:
            compact = line.upper().replace(' ', '')
            if 'EIGENVALUEOUTPUT' in compact:
                in_table = True
                continue
            if not in_table:
                continue
            items = line.split()
            if len(items) == 2:
                try:
                    mode = int(items[0])
                    value = float(items[1])
                except ValueError:
                    continue
                eigenvalues.append((mode, value))
            elif eigenvalues and not line.strip():
                # A blank line after the first entries ends the table.
                in_table = False

    return eigenvalues


def get_nodal_field(frame, name='U'):
    '''
    Get nodal field data as {instance name: (labels [n], data [n,ncomp])}.
    '''
    if name not in frame.fieldOutputs.keys():
        return {}

    field = frame.fieldOutputs[name]
    labels_of_instance = {}
    values_of_instance = {}

    for block in field.bulkDataBlocks:
        if block.nodeLabels is None or len(block.nodeLabels) == 0:
            continue
        name_instance = 'ASSEMBLY'
        if block.instance is not None:
            name_instance = str(block.instance.name)
        labels_of_instance.setdefault(name_instance, []).append(
            np.asarray(block.nodeLabels, dtype=np.int64))
        values_of_instance.setdefault(name_instance, []).append(
            np.asarray(block.data, dtype=np.float64))

    output = {}
    for name_instance in labels_of_instance:
        labels = np.concatenate(labels_of_instance[name_instance])
        data = np.concatenate(values_of_instance[name_instance], axis=0)
        output[name_instance] = (labels, data)

    return output


def copy_mesh(odb_source, odb_slim):
    '''
    Copy the mesh of every instance into a new output database.
    '''
    instances = {}

    for name, instance in odb_source.rootAssembly.instances.items():

        part = odb_slim.Part(name=str(name), embeddedSpace=THREE_D,
                             type=DEFORMABLE_BODY)

        node_labels = []
        node_coordinates = []
        for node in instance.nodes:
            node_labels.append(int(node.label))
            node_coordinates.append(tuple([float(x) for x in node.coordinates]))
        part.addNodes(labels=node_labels, coordinates=node_coordinates)

        element_labels = {}
        element_connectivity = {}
        for element in instance.elements:
            element_type = str(element.type)
            element_labels.setdefault(element_type, []).append(int(element.label))
            element_connectivity.setdefault(element_type, []).append(
                tuple([int(i) for i in element.connectivity]))
        for element_type in element_labels:
            part.addElements(labels=element_labels[element_type],
                             connectivity=element_connectivity[element_type],
                             type=element_type)

        instance_slim = odb_slim.rootAssembly.Instance(name=str(name), object=part)
        instances[str(name)] = instance_slim

        # The component sets (covers, spars, stringers, ribs) are what makes the
        # slim database navigable in the viewer, so keep them.
        n_sets = 0
        for name_set in instance.elementSets.keys():
            labels = [int(e.label) for e in instance.elementSets[name_set].elements]
            if not labels:
                continue
            instance_slim.ElementSetFromElementLabels(
                name=str(name_set), elementLabels=labels)
            n_sets += 1

        print('    [slim odb] %-10s nodes=%d elements=%d sets=%d'
              % (name, len(node_labels),
                 sum([len(v) for v in element_labels.values()]), n_sets))

    return instances


def add_nodal_field(frame, instances, data_of_instance, name, description):
    '''
    Add a nodal vector field to a frame of the slim output database.
    '''
    field = frame.FieldOutput(name=name, description=description, type=VECTOR,
                              componentLabels=('U1', 'U2', 'U3'),
                              validInvariants=(MAGNITUDE, ))
    for name_instance, (labels, data) in data_of_instance.items():
        if name_instance not in instances:
            continue
        field.addData(position=NODAL, instance=instances[name_instance],
                      labels=[int(i) for i in labels],
                      data=[tuple([float(x) for x in row]) for row in data])
    return field


def add_element_field(frame, instances, envelope, name, description):
    '''
    Add a scalar element field to a frame of the slim output database.
    '''
    field = frame.FieldOutput(name=name, description=description, type=SCALAR)
    for name_instance, (labels, values) in envelope.items():
        if name_instance not in instances:
            continue
        data = [(float(v), ) for v in values]
        try:
            field.addData(position=WHOLE_ELEMENT, instance=instances[name_instance],
                          labels=[int(i) for i in labels], data=data)
        except Exception:
            field.addData(position=CENTROID, instance=instances[name_instance],
                          labels=[int(i) for i in labels], data=data)
    return field


def main(argv):

    arguments = parse_arguments(argv)
    name_job = arguments['job']
    path_odb = arguments['odb']

    if not os.path.isfile(path_odb):
        raise IOError('Output database not found: %s' % path_odb)

    print('>>> --------------------')
    print('    [postprocess] %s' % path_odb)

    odb = openOdb(path=path_odb, readOnly=True)

    name_static = get_static_step(odb)
    name_buckle = get_buckle_step(odb)

    #* A pure buckling analysis has no loaded state, so there is no failure
    #* index to reduce and only the mode shape is kept.
    frame_static = None
    variables = []
    if name_static is not None:
        frame_static = odb.steps[name_static].frames[-1]
        variables = get_failure_variable_names(frame_static)
        print('    static step   = %s (frame %d)'
              % (name_static, frame_static.incrementNumber))
    else:
        print('    static step   = None (buckling only)')
    print('    buckling step = %s' % name_buckle)
    print('    failure index fields = %s' % (variables, ))

    #* ============================================
    #* Through-thickness envelope of the failure indices
    #* ============================================
    envelopes = {}
    summary = {
        'job': name_job,
        'static_step': name_static,
        'buckling_step': name_buckle,
        'variables': variables,
        'global_max': {},
        'max_of_instance': {},
    }

    #* The von Mises envelope of the stress field, if one was requested, so
    #* that metallic parts can be checked in the same way.
    name_stress = None
    if frame_static is not None:
        name_stress = get_stress_field_name(frame_static)
    if name_stress is not None:
        variables = list(variables) + ['S_MISES']
        summary['variables'] = variables
        print('    stress field         = %s (von Mises envelope)' % name_stress)

    arrays = {}
    for variable in variables:
        if variable == 'S_MISES':
            field = frame_static.fieldOutputs[name_stress].getScalarField(
                invariant=MISES)
        else:
            field = frame_static.fieldOutputs[variable]
        envelope = compute_element_envelope(field)
        envelopes[variable] = envelope

        global_max = 0.0
        max_of_instance = {}
        for name_instance, (labels, values) in envelope.items():
            arrays['%s__%s__labels' % (name_instance, variable)] = labels
            arrays['%s__%s' % (name_instance, variable)] = values
            if values.size == 0:
                continue
            index = int(np.argmax(values))
            max_of_instance[name_instance] = {
                'max': float(values[index]),
                'element': int(labels[index]),
                'n_elements': int(values.size),
            }
            global_max = max(global_max, float(values[index]))

        summary['global_max'][variable] = global_max
        summary['max_of_instance'][variable] = max_of_instance
        print('    max|%-9s| = %12.6g  (%d instances)'
              % (variable, global_max, len(envelope)))

    summary['eigenvalues'] = [
        {'mode': mode, 'eigenvalue': value}
        for mode, value in read_eigenvalues(name_job + '.dat')]
    if summary['eigenvalues']:
        print('    buckling eigenvalues = %s'
              % [round(item['eigenvalue'], 6) for item in summary['eigenvalues']])

    fname_npz = '%s_failure_envelope.npz' % name_job
    fname_json = '%s_failure_summary.json' % name_job
    np.savez_compressed(fname_npz, **arrays)
    with open(fname_json, 'w') as f:
        json.dump(summary, f, indent=2)
    print('    [write] %s (%.2f MB)'
          % (fname_npz, os.path.getsize(fname_npz)/1024.0**2))
    print('    [write] %s' % fname_json)

    #* ============================================
    #* Slim output database for visualisation
    #* ============================================
    if arguments['slim_odb']:

        path_slim = '%s_slim.odb' % name_job
        if os.path.isfile(path_slim):
            os.remove(path_slim)

        odb_slim = odbAccess.Odb(
            name=name_job, analysisTitle='%s (reduced)' % name_job,
            description='First buckling mode, deformed shape, and the '
                        'through-thickness envelope of the failure indices.',
            path=path_slim)

        instances = copy_mesh(odb, odb_slim)

        if name_buckle is not None and len(odb.steps[name_buckle].frames) > 1:
            frames = odb.steps[name_buckle].frames
            frame_mode = frames[1]
            step_slim = odb_slim.Step(name=str(name_buckle),
                                      description='Linear perturbation buckling',
                                      domain=TIME, timePeriod=1.0)
            frame_slim = step_slim.Frame(
                incrementNumber=1, frameValue=float(frame_mode.frameValue),
                description=str(frame_mode.description))
            add_nodal_field(frame_slim, instances,
                            get_nodal_field(frame_mode, 'U'),
                            'U', 'Buckling mode shape')
            print('    [slim odb] %s: %s' % (name_buckle, frame_mode.description))

        step_slim = odb_slim.Step(name=str(name_static),
                                  description='Static loading',
                                  domain=TIME, timePeriod=1.0)
        frame_slim = step_slim.Frame(
            incrementNumber=int(frame_static.incrementNumber),
            frameValue=float(frame_static.frameValue),
            description='Deformed shape and failure index envelope')

        add_nodal_field(frame_slim, instances,
                        get_nodal_field(frame_static, 'U'),
                        'U', 'Spatial displacement')

        for variable in variables:
            add_element_field(
                frame_slim, instances, envelopes[variable],
                '%s_MAX' % variable,
                'Through-thickness max|%s| of each element' % variable)

        odb_slim.save()
        odb_slim.close()
        print('    [write] %s (%.2f MB, source %.2f MB)'
              % (path_slim, os.path.getsize(path_slim)/1024.0**2,
                 os.path.getsize(path_odb)/1024.0**2))

    odb.close()

    if arguments['delete_source_odb']:
        if not arguments['slim_odb']:
            raise ValueError(
                'Refusing to delete %s without --slim-odb: the buckling mode '
                'and the deformed shape would be lost.' % path_odb)
        os.remove(path_odb)
        print('    [delete] %s' % path_odb)

    print('>>>')


if __name__ == '__main__':
    main(sys.argv)
