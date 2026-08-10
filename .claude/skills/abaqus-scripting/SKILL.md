---
name: abaqus-scripting
description: Write or modify Abaqus/Python scripts in this repository — building models with AbaqusTools (Model, Part, NodeOperation, OdbOperation, PBC/LBC, materials, LaRC05), looking up the Abaqus kernel API in the vendored abqpy stubs, editing *.inp files, and launching abaqus cae. Use whenever the task touches AbaqusTools/, wingbox/, example/, a *.inp/*.odb file, or any Abaqus scripting question.
---

# Abaqus scripting in this repository

## 1. The environment, before anything else

| | |
|---|---|
| Solver | **Abaqus 2023** |
| Kernel interpreter | **Python 2.7** — everything under `AbaqusTools/`, `wingbox/` and `example/` that runs inside Abaqus must be Python 2.7 compatible |
| Launch | `abaqus cae noGUI=script.py` (batch) or `abaqus cae script=script.py` (with GUI) |
| Script arguments | after a bare `--`, e.g. `abaqus cae noGUI=wingbox_model.py -- --params my-parameters.json`; parsed by hand in `wingbox/params.py` because `argparse` is awkward under 2.7 |
| Units | **N-mm-tonne** throughout: mm, N, MPa, tonne/mm³, fracture energy in N/mm (= kJ/m²) |

**Python 2.7 rules.** No f-strings, no `pathlib`, no `subprocess.run`, no type annotations in
signatures, no dict/set comprehension features beyond 2.7, `class Foo(object)` not `class Foo`,
`print('%s'%(x))` not `print(f'{x}')`, integer division is `//` when you mean it. The type checker
cannot enforce this (pyright only understands Python 3) — it is on you.
Post-processing scripts that run *outside* Abaqus (numpy/scipy/matplotlib) may use Python 3.

**The `IS_ABAQUS` guard.** `AbaqusTools/__init__.py` sniffs `sys.argv` for `-tmpdir` / `-cae` /
`-interactive` and only then imports the Abaqus modules. Every module that touches the Abaqus API
must repeat the pattern, or it becomes unimportable outside Abaqus:

```python
from AbaqusTools import IS_ABAQUS

if IS_ABAQUS:
    from abaqus import *
    from abaqusConstants import *
    from caeModules import *
```

Do not remove this guard. See §2 for why it matters even more now that abqpy is vendored.

## 2. abqpy — read it, do not import it

`external/abqpy` is a git submodule pinned to abqpy's **2023** branch, matching the solver version.
It is 1400+ `.py` files carrying the full Abaqus kernel API with the official documentation as
docstrings — `external/abqpy/src/abaqus/**`, plus `external/abqpy/src/abaqusConstants.py`.

**Use it as the API reference.** When you are unsure of a keyword argument, a symbolic constant or a
return type, *read the stub* instead of guessing:

```bash
grep -rn "def StaticStep" external/abqpy/src/abaqus/Step/
grep -rn "def HashinDamageInitiation" external/abqpy/src/abaqus/Material/
sed -n '/class Odb/,/def /p' external/abqpy/src/abaqus/Odb/Odb.py
```

This is the single highest-value habit in this repo: it turns "what does `ContactStd` accept" from a
recall problem into a lookup.

**Never `import abaqus` outside Abaqus.** abqpy is not a passive stub package. `abaqus/__init__.py`
calls `abqpy.run(cae=True)` at import time, which tries to **re-launch your script through the real
`abaqus cae` executable**. Importing `abaqusConstants` triggers it too (it re-exports from
`abaqus.UtilityAndView`). Verified behaviour:

```
>>> import abaqusConstants
RuntimeError: Cannot find the main script file, please run the script in a file.
```

So the `IS_ABAQUS` guard is not a nicety, it is what stops a stray `python postprocess.py` from
spawning CAE. Keep the stubs on the *type checker's* path (`pyrightconfig.json` → `extraPaths`),
not on the runtime path. Installing abqpy with pip is unnecessary for this repo and is only worth it
if you specifically want abqpy's own CLI runner.

Setup instructions live in `external/README.md`.

## 3. What to reuse from this repository

Reach for these before writing raw Abaqus API calls. Everything below already exists and is used by
the scripts in `example/` and `wingbox/`.

**`AbaqusTools/model.py` — `Model`**
Subclass it and fill in the hooks; `Model.build()` calls them in this fixed order:
`setup_parts` → `setup_assembly` → `setup_steps` → `setup_outputs` → `setup_interactions` →
`setup_loads` → `setup_jobs`. `setup_property` runs earlier, from `initialization()`.
Useful members: `create_geometry_set`, `create_node_set`, `get_assembly_set`, `get_nodes_in_set`,
`create_reference_point` / `get_reference_point` / `create_reference_point_set`,
`create_analytical_mapped_field`, `create_contact_constraints`, `submit_job`, `write_job_inp`,
`save_cae`, `set_view`.

**`AbaqusTools/model.py` — `NodeOperation`** (static methods, also the base of the BC classes)
`get_set`, `get_nodes_from_face`, `exclude_forbidden_nodes`, `create_node_set`,
`create_face_node_set`, `create_paired_node_sets`. This is where node sorting and pairing for
boundary conditions lives — do not re-derive it.

**`AbaqusTools/part.py` — `Part`**
Subclass and fill in; `Part.build()` calls `create_sketch` → `create_part` → `create_partition` →
`create_surface` → `create_set` → (`set_seeding` → `create_mesh` → `set_element_type` →
`set_section_assignment` → `set_composite_layups`, skipped when `is_only_geometry`).
Geometry pickers: `get_vertex/edge/face/cell` (one object) and `get_vertices/edges/faces/cells`
(arrays). Also `create_datum_point`, `create_datum_csys_3p`, `get_vertices_on_face`,
`get_CompositeLayup_thickness`.

**`AbaqusTools/materials.py`** — material data as plain dicts; see §4.

**`AbaqusTools/pbc.py`, `AbaqusTools/lin_bc.py`** — `PBC_Beam`, `PBC_3DOrthotropic`,
`LBC_3DOrthotropic`, `LBC_InPlaneLoad`. Constraint-equation generation, node pairing and
engineering-constant back-out for homogenisation. Substantial, validated code — extend, never
reimplement.

**`AbaqusTools/odb.py` — `OdbOperation`** — ODB reading, including the `fieldOutput` index ↔ element
/node label mapping (`convert_IdxFO_to_Label`), which is the part of the ODB API most likely to be
silently wrong if hand-rolled.

**`AbaqusTools/larc05.py`** — LaRC05 failure criterion in Python (`UVARM`, `PlyProperty`,
`FailureCriteria`), mirroring the Fortran subroutines in `LaRC05/`.

**`AbaqusTools/functions.py`** — `LayupParameters.candidate_composite_layup` (stacking-sequence
design), `load_parameters`, `clean_pyc_files`, `clean_temporary_files`.

## 4. Materials are data, not methods

Material cards live in `AbaqusTools/materials.py::MATERIAL_LIBRARY` (`Steel`, `Ti-6Al-4V`,
`Aluminum-7075`, `IM7/8551-7`), stored in the N-mm-tonne system with a `stress_columns` entry per
table so `Model.create_material(..., unit_length='m')` can convert.

```python
model.create_material('Ti-6Al-4V')                                  # default elasticity
model.create_material('IM7/8551-7', elastic_type='ENGINEERING_CONSTANTS')
model.create_section('Ti-6Al-4V')                                   # HomogeneousSolidSection
```

**To add a material, add a dict entry — do not add a method.** The module imports nothing from
Abaqus on purpose, so it can be read and tested in a normal Python.

`create_material_IM785517()` is kept because it carries *policy*, not data: it reads
`pMesh['failure_model']` and attaches Hashin damage (`create_hashin_damage`) or announces the LaRC05
UMAT/UVARM path. `create_material_steel/titanium` and `create_section_*` are thin aliases kept so the
frozen scripts under `example/` keep running; new code should call `create_material`/`create_section`.

## 5. Conventions worth not breaking

- **Step name is `'Loading'`, previous is `'Initial'`.** `create_static_step()` pins this along with
  `minInc=1E-15`. `wingbox/` renames to `'Preload'`/`'Buckling'` for `analysis_type='static_buckle'`
  via `_get_static_step_name()`; follow that pattern rather than inventing new names.
- **Model name is `'Model-1'`** (`Model.name_model`) — Abaqus's default, deliberately not changed.
- **Reference points**: always use `create_reference_point(x, y, z, name_rp)`. It creates the feature
  and immediately renames `'RP-1'`, and raises if the name already exists. Abaqus's auto-generated
  `'RP-1'` must never be relied on.
- **Parameters** flow as three dicts, `pGeo` / `pMesh` / `pRun`, loaded from JSON. New knobs go into
  the JSON files, read with `.get(key, default)` so old parameter files keep working.
- **`temp*` is gitignored** and holds throwaway run directories with vendored copies of
  `AbaqusTools/`. Never edit those copies — they are stale by design. Edit `AbaqusTools/`.

## 6. `findAt` traps

The single most common source of silent breakage in this repo.

```python
myPrt.vertices.findAt((pt))     # -> Vertex object
myPrt.vertices.findAt((pt,))    # -> VertexArray sequence
```

One extra comma changes the return type, and most Abaqus APIs (`Set`, `Surface`, section
assignment) demand the **array** form. `Part.get_*` wraps this: pass `toArray=True` when the result
is going into a region argument.

`findAt` only matches within 1E-6 of the geometry. When a point is computed rather than exact, use
`getClosest=True` with an explicit `searchTolerance`:

```python
faces = Part.get_faces(myPrt, [(x, y, z)], getClosest=True, searchTolerance=1.0)
```

`findAt` returns objects whose indices are invalidated by any later partitioning or geometry edit.
Pick geometry *after* the part is fully partitioned, or re-pick.

## 7. Recipes

Snippets that used to be one-line wrapper methods. Inline them where needed instead of adding a
method to `Model`.

**Amplitudes**

```python
self.model.TabularAmplitude(name='Constant-Amp', timeSpan=STEP,
        smooth=SOLVER_DEFAULT, data=((0.0, 1.0), (1.0, 1.0)))
self.model.TabularAmplitude(name='Ramp-Amp', timeSpan=STEP,
        smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
```

**Explicit dynamic step** — element types must then come from the `EXPLICIT` `elemLibrary`.

```python
self.model.ExplicitDynamicsStep(
        name='Loading', previous='Initial',
        description='Dynamic (explicit) simulation',
        nlgeom=ON, improvedDtMethod=ON)
```

**Buckling step** — already implemented properly in `wingbox/wingbox_model.py::setup_steps`, which
assembles `BuckleStep` options (`numEigen`, `vectors`, `maxIterations`, `minEigen`, `blockSize`,
`maxBlocks`) from `pRun` and handles the `static_buckle` two-step case. Copy that, do not re-derive.

**Fixed camera for reproducible screenshots** — record the numbers from a CAE session you like
(`session.views` after orienting by hand), then hard-code them in the *model script*, not in
`AbaqusTools`:

```python
session.View(name='User-1', nearPlane=2549.1, farPlane=4259.2,
        width=1669.3, height=809, projection=PERSPECTIVE,
        cameraPosition=(-1797.7, 1551.9, 2152.4),
        cameraUpVector=(0.42948, 0.88149, -0.19628),
        cameraTarget=(771.82, -44.15, 607.24),
        viewOffsetX=-3.0688, viewOffsetY=-3.6258, autoFit=OFF)
session.viewports['Viewport: 1'].view.setValues(session.views['User-1'])
```

**Contact with friction** — the choices here are deliberate: penalty friction µ from
`pRun['contact_friction_coef']`, hard normal contact with separation allowed. See
`Model.create_interaction_property_contact` and `create_contact_constraints`; the secondary surface
of a tie pair must be the smaller one.

## 8. Editing `*.inp` files

Some things have no kernel API and must be done by rewriting the input deck after
`write_job_inp()` and before submitting: `Model.write_static_step_inp`,
`write_output_field_frequency_interval`, `write_IM785517_property_table_inp`.

All of them scan the deck line by line, match a keyword, write replacement lines, and drop the
original block. **How the original block is dropped is the part that goes wrong**, and there are two
correct answers:

*Fixed line count* — valid only when the replaced block has a known, constant length.
`write_static_step_inp` skips `N_SKIP_LINE_STATIC = 2` after `*Static` (the keyword plus its single
data line); `write_output_field_frequency_interval` skips 1 after `*Output, field`.

*Scan to the end of the block* — required when the length varies.
`write_IM785517_property_table_inp` must do this: after `*Material, name=IM7/8551-7` it keeps
skipping while the line is a material option (`*ELASTIC`, `*DENSITY`) or one of its data lines, and
stops at the first keyword that is not. A fixed count of 3 was wrong here, because
`*ELASTIC, type=ENGINEERING CONSTANTS` spans two data lines; the leftover line made Abaqus abort
with `TOO FEW LINES TO DEFINE THIS MATERIAL OPTION`. Do not reintroduce a fixed count for any block
whose data lines depend on the material definition.

Other details worth keeping: keyword matching is on `line.split()`, so it tolerates whitespace but
is **case-sensitive** — Abaqus writes `*Material, name=IM7/8551-7` with exactly that capitalisation,
while the block-end test upper-cases the keyword before comparing. Always print what was overwritten;
these edits are otherwise invisible until the solver complains.

## 9. Running LaRC05 jobs

The LaRC05 path (`pMesh['failure_model'] == 'larc05'`) does not submit from inside CAE. The deck is
written, patched (§8), and only then submitted from `wingbox/run.py`:

```bash
abaqus interactive job=<name> user=uvarm.f90 cpus=<n> standard_parallel=solver
```

- **`standard_parallel=solver` is mandatory.** The LaRC05 user subroutine keeps module-level state
  and is not thread-safe; a parallel element loop segfaults. This flag keeps the element loop serial
  while the linear solver still uses every cpu. Never drop it to "speed things up".
- **`abaqus_v6.env` must survive cleanup.** `clean_temporary_files()` deletes job scratch `*.env`
  files but explicitly spares `abaqus_v6.env`, which carries the local Fortran compiler override
  between the model build and the job submission. Deleting it makes the user subroutine fail to
  compile.
- The Fortran sources live in `LaRC05/` and are tracked deliberately, since the run cannot be
  reproduced without them.

## 10. Before finishing

- Python 2.7 syntax in everything that runs inside Abaqus.
- Guarded Abaqus imports (`if IS_ABAQUS:`).
- API kwargs checked against `external/abqpy/src/`, not from memory.
- Arrays where regions are expected (`toArray=True`).
- New material → `MATERIAL_LIBRARY` entry; new tunable → JSON + `.get(key, default)`.
- Changes go in `AbaqusTools/`, never in a `temp*/AbaqusTools/` copy.
- You cannot run Abaqus yourself to verify. Say plainly what was checked (syntax, imports, data) and
  what still needs a real solver run.
