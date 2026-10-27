# External dependencies

## `abqpy` — Abaqus/Python type hints

[abqpy](https://github.com/haiiliin/abqpy) provides the whole Abaqus kernel API as annotated Python
modules, with the official documentation as docstrings. It is vendored here as a git submodule so
that the API reference is available offline, greppable, and pinned to a known version.

**The submodule tracks abqpy's `2023` branch, matching Abaqus 2023.** abqpy is versioned by Abaqus
release; a mismatched branch will suggest keyword arguments and Python 3 syntax that the Abaqus 2023
kernel (Python 2.7) rejects. If you upgrade Abaqus, switch the branch in `.gitmodules` and in the
submodule, then commit the new pointer.

### Setup

Fresh clone of this repository:

```bash
git clone --recurse-submodules --shallow-submodules <this-repo>
```

Already cloned without submodules:

```bash
git submodule update --init --depth 1 external/abqpy
```

That is the whole setup. Nothing to install, nothing to add to `requirements.txt`.

### How the type hints reach the editor

`pyrightconfig.json` at the repository root adds `external/abqpy/src` to `extraPaths`, so
Pylance/pyright resolve `abaqus`, `abaqusConstants`, `caeModules`, `odbAccess` and give completions
and signatures for them. `typeCheckingMode` is `off`: we want the hints, not thousands of
diagnostics on Python-2-era code.

VS Code picks this up automatically. Reload the window after the first `git submodule update`.

### Do not import abqpy at runtime

abqpy is **not** a passive stub package. `external/abqpy/src/abaqus/__init__.py` calls
`abqpy.run(cae=True)` on import, which tries to re-launch the current script through the real
`abaqus cae` executable. `abaqusConstants` re-exports from `abaqus.UtilityAndView`, so it triggers
the same thing:

```
$ python -c "import abaqusConstants"
RuntimeError: Cannot find the main script file, please run the script in a file.
```

This is why `AbaqusTools/__init__.py` guards every Abaqus import behind `IS_ABAQUS`, and why
`external/abqpy/src` is on the type checker's path and not on `sys.path`. Keep it that way.

There is no reason to `pip install abqpy` for this repository. If you want abqpy's own CLI runner
for some other project, install it in a separate environment:

```bash
pip install abqpy==2023.*
```

### Updating the pin

```bash
git -C external/abqpy fetch --depth 1 origin 2023
git -C external/abqpy checkout FETCH_HEAD
git add external/abqpy && git commit -m "Bump abqpy pin"
```
