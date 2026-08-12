# -*- coding: utf-8 -*-
'''
Fast extraction of the wingbox results an optimizer needs.

Unlike `postprocess_failure.py`, this script does NOT scan the section
point fields of the output database (minutes for a 100-ply model).  It
collects the three quantities of a design evaluation from cheap sources:

- maximum LaRC05 failure index of the whole structure — read from
  `larc05_fi_track_<pid>.txt`, the per-rank tracking files rewritten by
  UVARM in `uvarm.f90` (module larc05Track) whenever a rank's running
  maximum failure index increases (a legacy single `larc05_fi_track.txt`
  is also accepted),
- buckling eigenvalues — parsed from the `<job>.dat` text file,
- wing tip maximum displacement — the nodal U field of the last static
  frame, the only field read from the output database.

The summary is written to `<job>_failure_summary.json` with a schema
compatible with `postprocess_failure.py` (the fast mode adds
`max_fi_location` and omits `max_of_instance`), so downstream scripts can
consume either one.

Usage (Abaqus Python, no CAE license needed)
--------------------------------------------
    abaqus python extract_results.py --job Job_WB
'''
import os
import sys
import json

import numpy as np

import odbAccess
from odbAccess import openOdb

from postprocess_failure import (parse_arguments, get_static_step,
                                 get_buckle_step, read_eigenvalues,
                                 compute_tip_displacement)

#* Written by UEXTERNALDB in uvarm.f90, one file per MPI rank
#* (larc05_fi_track_<pid>.txt), one record per increment:
#* step, increment, running rank-local max FI, element, layer.
TRACK_GLOB = 'larc05_fi_track_*.txt'
TRACK_FILE_LEGACY = 'larc05_fi_track.txt'


def read_tracked_max_fi(pattern=TRACK_GLOB):
    '''
    Merge the UVARM tracking files of all MPI ranks.

    Each rank accumulates the maximum failure index of its own elements,
    so the global maximum is the largest last-record value over all
    files.  A legacy single `larc05_fi_track.txt` is also accepted.

    Returns
    ----------------
    result: dict, None
        `max_fi` (uncapped global maximum failure index), `element`,
        `layer`, `step` and `increment` of the winning record; None when
        no tracking file is found.
    '''
    import glob

    fnames = sorted(glob.glob(pattern))
    if os.path.isfile(TRACK_FILE_LEGACY):
        fnames.append(TRACK_FILE_LEGACY)

    best = None
    for fname in fnames:
        with open(fname, 'r') as f:
            for line in f:
                items = line.split()
                if len(items) < 3:
                    continue
                try:
                    record = {
                        'step': int(items[0]),
                        'increment': int(items[1]),
                        'max_fi': float(items[2]),
                        'element': int(items[3]) if len(items) > 3 else 0,
                        'layer': int(items[4]) if len(items) > 4 else 0,
                    }
                except ValueError:
                    continue
                if best is None or record['max_fi'] > best['max_fi']:
                    best = record

    return best


def compute_tip_displacement_history(odb, name_static):
    '''
    Maximum displacement magnitude of the wing tip from the history output.

    The model writes U of the tip cross-section nodes as history output
    (`Tip-Output` request, one history region per node, named
    'Node <INSTANCE>.<LABEL>'), so only the last frame of a few hundred
    nodes is read instead of the full nodal field.

    Returns the same dict as `compute_tip_displacement` (without `z_tip`,
    which the history data does not carry), or None when no tip history
    regions are present.
    '''
    if name_static is None:
        return None

    result = None
    for name_region, region in odb.steps[name_static].historyRegions.items():
        outputs = region.historyOutputs
        if 'U1' not in outputs.keys():
            continue
        u1 = outputs['U1'].data
        u2 = outputs['U2'].data
        u3 = outputs['U3'].data
        if len(u1) == 0:
            continue
        U = [float(u1[-1][1]), float(u2[-1][1]), float(u3[-1][1])]
        magnitude = float(np.sqrt(u1[-1][1]**2 + u2[-1][1]**2 + u3[-1][1]**2))

        # 'Node LOFTING.222' -> instance LOFTING, node 222
        tokens = str(name_region).split()
        instance, node = tokens[-1].rsplit('.', 1)
        candidate = {
            'max_magnitude': magnitude,
            'node': int(node),
            'instance': instance,
            'z_tip': None,
            'U': U,
        }
        if result is None or magnitude > result['max_magnitude']:
            result = candidate

    return result


def main(argv):

    arguments = parse_arguments(argv)
    name_job = arguments['job']
    path_odb = arguments['odb']

    print('>>> --------------------')
    print('    [extract] %s' % name_job)

    #* Maximum failure index from the UVARM tracking file.
    tracked = read_tracked_max_fi()
    if tracked is None:
        raise IOError('UVARM tracking file not found or empty: %s — was '
                      'the job run with the runfile uvarm.f90?' % TRACK_GLOB)
    print('    max FI (UVARM6, uncapped) = %.6g (element %d, layer %d)'
          % (tracked['max_fi'], tracked['element'], tracked['layer']))

    #* Buckling eigenvalues from the data file.
    eigenvalues = [{'mode': mode, 'eigenvalue': value}
                   for mode, value in read_eigenvalues(name_job + '.dat')]
    if eigenvalues:
        print('    buckling eigenvalues = %s'
              % [round(item['eigenvalue'], 6) for item in eigenvalues])

    #* Wing tip displacement: history output of the tip node set when
    #* available (cheap), otherwise the full nodal U field of the last
    #* static frame (slow fallback for decks without the tip request).
    tip = None
    name_static = None
    name_buckle = None
    if os.path.isfile(path_odb):
        odb = openOdb(path=path_odb, readOnly=True)
        name_static = get_static_step(odb)
        name_buckle = get_buckle_step(odb)
        tip = compute_tip_displacement_history(odb, name_static)
        if tip is not None:
            print('    tip displacement   = %.6g mm (node %d of %s, history)'
                  % (tip['max_magnitude'], tip['node'], tip['instance']))
        else:
            frame_static = (odb.steps[name_static].frames[-1]
                            if name_static is not None else None)
            tip = compute_tip_displacement(odb, frame_static)
            if tip is not None:
                print('    tip displacement   = %.6g mm (node %d of %s, z = %.1f mm)'
                      % (tip['max_magnitude'], tip['node'], tip['instance'],
                         tip['z_tip']))
        odb.close()
    else:
        print('    [warn] output database not found: %s' % path_odb)

    summary = {
        'job': name_job,
        'static_step': name_static,
        'buckling_step': name_buckle,
        'variables': ['UVARM6'],
        'global_max': {'UVARM6': tracked['max_fi']},
        'max_fi_location': {
            'element': tracked['element'],
            'layer': tracked['layer'],
        },
        'eigenvalues': eigenvalues,
        'tip_displacement': tip,
    }

    fname_json = '%s_failure_summary.json' % name_job
    with open(fname_json, 'w') as f:
        json.dump(summary, f, indent=2)
    print('    [write] %s' % fname_json)
    print('>>>')


if __name__ == '__main__':
    main(sys.argv)
