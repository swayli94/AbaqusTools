# Wing box with S4R elements

## Running a case

`run.py` builds the model, runs the job, and reduces the output database.
It never reads the parameter file by a fixed name any more:

```bash
python -u run.py                                  # parameters-simple-wingbox.json
python -u run.py --params my-parameters.json      # any other file
./run.sh parameters-inner-wingbox.json            # same, after cleaning the previous run
```

`run.sh` copies nothing: it expects the working folder to already hold
`AbaqusTools/`, the python files of this folder, the airfoil `*.dat`, the
parameter file, and — only for LaRC05 — the `*.f90` of `LaRC05/`.  The two
parameter files shipped here are `parameters-simple-wingbox.json` (a small
three-section box that runs in minutes) and `parameters-inner-wingbox.json`
(the full inner-wing case).

`run.py` forwards the file to the model script, so calling Abaqus directly is
equivalent:

```bash
abaqus cae noGUI=wingbox_model.py -- --params my-parameters.json
```

Abaqus/CAE submits jobs asynchronously.  `"wait_for_completion"` therefore
defaults to `true`: the CAE session waits for the solver, and only then does
`run.py` remove the temporary files (`*.com`, `*.env`, `*.sim`, `*.stt`), which
the solver is still using while it runs.  Setting it to `false` keeps those
files in place, but then nothing waits for the job.

## What is analysed

| parameter | values | meaning |
|---|---|---|
| `pRun["analysis_type"]` | `"static"`, `"buckle"`, `"static_buckle"` | whether a buckling step is solved, and whether a general static step is solved |
| `pMesh["failure_model"]` | `"Hashin"`, `"LaRC05"`, `"none"` / absent / `null` | which failure criterion is defined and written; `"none"` switches failure analysis off entirely |
| `pRun["output_rib_variables"]` | list of variables, `[]` to disable | stress output of the ribs, defaulting to `["S"]` for metallic ribs |

A metallic rib (`pMesh["rib_material_type"]: "aluminum"`) has no composite
layup, so it appears in no layup output request.  Its strength can only be
checked from its stress field, which is why `Rib-Output-rib_*` requests are
created by default in that case, and `postprocess_failure.py` reduces them to a
per-element `S_MISES_MAX`.

## Output database size

The failure indices are evaluated at every section point of every ply, so a
100-ply cover writes 300 values per element with a three-point Simpson rule.
Two settings control how much of that reaches disk.

`pRun["layup_output_ply_locations"]` selects the section points written by the
layup output requests:

| value                | section points per ply | note                                    |
|----------------------|------------------------|-----------------------------------------|
| `"all"` (default)    | all (3 for Simpson)    | rigorous                                |
| `"top_bottom"`       | 2                      | keeps the ply extremes, ~30% smaller    |
| `"mid"`              | 1                      | smallest, least conservative            |

`postprocess_failure.py` then reduces the raw database to what a failure
analysis needs.  For every element it takes `max(|value|)` over all section
points of each failure index, i.e., the through-thickness envelope, and writes

- `<job>_failure_envelope.npz` — element labels and envelope per instance,
- `<job>_failure_summary.json` — global/per-instance maxima, the buckling
  eigenvalues read from `<job>.dat`, and `tip_displacement`, the maximum
  displacement magnitude of the wing tip nodes (largest z coordinate) in the
  loaded state,
- `<job>_slim.odb` — the mesh, the first buckling mode shape, the deformed
  shape of the static step, and the envelope as `<VARIABLE>_MAX` element
  fields.

The slim database is typically one to two orders of magnitude smaller than the
source, which can then be deleted:

```json
"postprocess": {
    "enabled": true,
    "slim_odb": true,
    "delete_source_odb": false
}
```

`delete_source_odb` requires `slim_odb`, so the mode shape and the deformed
shape are never lost.  The script can also be run on its own:

```bash
abaqus python postprocess_failure.py --job Job_WB --slim-odb
```

Both the built-in Hashin criteria (`HSNFTCRT`, `HSNFCCRT`, `HSNMTCRT`,
`HSNMCCRT`) and the LaRC05 user-defined output variables (`UVARM1`, `UVARM2`,
...) are recognised, so the same post-processing applies after switching
`pMesh["failure_model"]` to `LaRC05`.
