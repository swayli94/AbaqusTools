# AbaqusTools

A comprehensive Python package for streamlining and automating ABAQUS finite element analysis workflows. This library provides high-level wrapper functions and object-oriented interfaces that make ABAQUS Python scripting more intuitive, efficient, and maintainable.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

ABAQUS is a powerful finite element analysis software widely used in engineering for structural, thermal, and multiphysics simulations. While ABAQUS provides extensive Python scripting capabilities, writing complex automation scripts often requires deep knowledge of its API structure and involves repetitive coding patterns.

AbaqusTools addresses these challenges by providing intuitive Python classes and functions that wrap the complex ABAQUS API, making it easier to:

- Create and manage ABAQUS models programmatically
- Automate repetitive modeling tasks
- Process simulation results efficiently
- Work with composite materials and layup designs
- Apply periodic boundary conditions (PBC) and linear boundary conditions (LBC)
- Evaluate composite failure with the LaRC05 criterion
- Run batch analyses and parameter studies

## Key Features

### Model Management

- Simplified interfaces for creating, modifying, and managing ABAQUS models
- Object-oriented `Model` class that handles Assembly, Step, Interaction, Load, and Job operations
- Intelligent environment detection that works both within and outside ABAQUS

### Part Operations

- Advanced part creation and manipulation with the `Part` class
- Geometric operations and automated mesh generation
- Operations in Sketch, Part, Property, and Mesh modules

### Composite Material Support

- Specialized `LayupParameters` class for composite layup design
- Automated stacking sequence generation following industry best practices
- Support for balanced, symmetric, and quasi-isotropic layups
- Built-in mechanical joint design guidelines

### Results Processing

- Streamlined `OdbOperation` class for ODB (Output Database) operations
- Extract and process simulation results efficiently
- Node and element data manipulation tools

### Periodic Boundary Conditions

- `PeriodicBC` base class with specialized subclasses: `PBC_Beam` (beam homogenization) and `PBC_3DOrthotropic` (3D RVE homogenization)
- Automated constraint equation generation using ABAQUS `*EQUATION` keyword
- Node pairing with automatic exclusion of forbidden nodes (corners/edges)

### Linear Boundary Conditions

- `LBC_3DOrthotropic` class for displacement-driven homogenization of 3D orthotropic RVEs
- Applies linear deformation boundary conditions as an alternative to PBC
- Automated node-set creation on opposite faces for engineering-constant extraction

### Composite Failure Criteria

- `UVARM` class implementing the LaRC05 failure criterion in Python (translated from Fortran subroutine)
- `PlyProperty` for storing fiber and matrix material properties (e.g. IM7/8551-7)
- `FailureCriteria` for matrix cracking, splitting, fiber tension, and fiber kinking indices
- Compatible with both 2D and 3D element formulations

### Utility Functions

- Parameter loading and management
- File cleanup utilities for temporary and compiled files
- Cross-platform compatibility tools

## Installation

### Prerequisites

For the best development experience, install `abqpy` alongside AbaqusTools to get comprehensive type hints for ABAQUS Python scripting:

```bash
pip install abqpy
```

This enables you to write ABAQUS Python scripts with full IDE support, code completion, and type checking, even without opening ABAQUS/CAE.

### Installing AbaqusTools

Clone this repository and install:

```bash
git clone https://github.com/your-username/AbaqusTools.git
cd AbaqusTools
# Add to your Python path or install as needed
```

## Quick Start

### Basic Model Creation

```python
from AbaqusTools import Model, Part

# Create a new model
model = Model(pGeo=geometry_params, pMesh=mesh_params, pRun=run_params)

# Create and build parts
part = Part(model, pGeo=geometry_params, pMesh=mesh_params)
part.build()

# Assemble, apply loads, and submit job
model.assembly()
model.step()
model.load()
model.job()
```

### Composite Layup Design

```python
from AbaqusTools.functions import LayupParameters

# Create layup parameters
layup = LayupParameters(THK_PLY=0.25)

# Generate candidate layups
candidate_layup = layup.candidate_composite_layup(n_ply=16, index=0)
```

### Extracting ODB Results

```python
from AbaqusTools import OdbOperation

# Open and process results
odb = OdbOperation('job_name')
# Extract field data, nodal results, etc.
```

### Applying Boundary Conditions

```python
from AbaqusTools.pbc import PBC_3DOrthotropic
from AbaqusTools.lin_bc import LBC_3DOrthotropic

# Periodic boundary conditions for 3D RVE homogenization
pbc = PBC_3DOrthotropic(model, instance)

# Or use linear boundary conditions
lbc = LBC_3DOrthotropic(model, instance)
```

## Running Scripts

When running AbaqusTools scripts with ABAQUS, use the following command structure:

```bash
# Clean temporary files and run script
abaqus cae script=your_script.py

# For Windows users with cleanup
clean.bat & abaqus cae script=your_script.py
```

> **Note**: Place your scripts in the same directory as the 'AbaqusTools' folder. It's recommended to clean `*.pyc` files before calling ABAQUS to avoid issues with cached bytecode.

## Documentation

Comprehensive documentation is available in the `docs/` directory, including:

- **Introduction and Getting Started**: Overview of AbaqusTools capabilities
- **Operation Guides**: Detailed guides for each ABAQUS module
- **API Reference**: Complete class and function documentation
- **Examples**: Sample scripts and usage patterns

To view the documentation locally:

```bash
cd docs/build/html
# Open index.html in your browser
```

## Project Structure

```text
AbaqusTools/
├── AbaqusTools/          # Main package
│   ├── __init__.py       # Package initialization and environment detection
│   ├── model.py          # Model and NodeOperation classes
│   ├── part.py           # Part class for geometry and meshing
│   ├── odb.py            # OdbOperation for results processing
│   ├── pbc.py            # PeriodicBC, PBC_Beam, PBC_3DOrthotropic classes
│   ├── lin_bc.py         # LBC_3DOrthotropic classes for linear boundary conditions
│   ├── functions.py      # Utility functions and LayupParameters
│   └── larc05.py         # LaRC05 failure criterion (UVARM, PlyProperty, FailureCriteria)
├── docs/                 # Sphinx documentation
├── example/              # Example scripts
└── LICENSE               # MIT License
```

## Environment Detection

AbaqusTools automatically detects whether it's running within an ABAQUS environment by checking command-line arguments. This allows the same code to work both in ABAQUS/CAE and in external Python environments for development and testing.

## Who Should Use AbaqusTools?

- **Simulation Engineers** automating repetitive modeling tasks
- **Research Scientists** conducting parameter studies and optimization
- **Composite Engineers** working with laminated materials and complex layups
- **CAE Analysts** looking to streamline ABAQUS workflows
- **Python Developers** building custom simulation tools and applications

## Dependencies

- Python (compatible with ABAQUS Python environment)
- NumPy
- SciPy
- abqpy (recommended for type hints)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## Author

**Runze LI** - Initial work and development

---

**Note**: This package is designed to work with ABAQUS finite element software. ABAQUS is a product of Dassault Systèmes SIMULIA Corp.
