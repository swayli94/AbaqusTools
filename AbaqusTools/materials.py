'''
Material data library.

This module holds **data only**: the material cards used by the models in this
repository, expressed as plain Python dictionaries. It deliberately imports
nothing from Abaqus, so it can be read, edited and unit-tested outside
Abaqus/CAE (e.g. by post-processing scripts running in a normal Python).

The Abaqus API calls that turn these dictionaries into `Material` and `Section`
objects live in `AbaqusTools.model.Model.create_material()` and
`Model.create_section()`.

Units
--------------
All tables are stored in the **N-mm-tonne** system, which is what every model in
this repository uses:

    length      mm
    force       N
    stress      MPa = N/mm^2
    density     tonne/mm^3
    energy      mJ = N*mm  (fracture energy: N/mm = kJ/m^2)

`Model.create_material(..., unit_length='m')` converts to the N-m-kg system on
the fly using `scale_table()`. The conversion needs to know which columns of a
table are stresses, hence the `stress_columns` entry next to every table.

Adding a material
--------------
Add an entry to `MATERIAL_LIBRARY`. No new method on `Model` is needed:

    >>> model.create_material('Ti-6Al-4V')
    >>> model.create_section('Ti-6Al-4V')
'''

#* Unit conversion factors from the N-mm-tonne system to the N-m-kg system.
UNIT_SCALE = {
    'mm': {'stress': 1.0, 'density': 1.0},
    'm':  {'stress': 1.0E6, 'density': 1.0E12},
}


#* Column layout of the Abaqus `Elastic` table for each elasticity type.
#* Used only as documentation for the `stress_columns` entries below.
ELASTIC_TABLE_COLUMNS = {
    'ISOTROPIC': ('E', 'Nu'),
    'LAMINA': ('E1', 'E2', 'Nu12', 'G12', 'G13', 'G23'),
    'ENGINEERING_CONSTANTS': ('E1', 'E2', 'E3', 'Nu12', 'Nu13', 'Nu23',
                              'G12', 'G13', 'G23'),
}


MATERIAL_LIBRARY = {

    'Steel': {
        'description': 'Generic structural steel with isotropic hardening',
        'density': 7.8E-9,
        'default_elastic': 'ISOTROPIC',
        'elastic': {
            'ISOTROPIC': {
                'type': 'ISOTROPIC',
                'table': ((2.1E5, 0.3), ),
                'stress_columns': (0, ),
            },
        },
        'plastic': {
            # (Yield stress, Plastic strain)
            'table': ((3.00E2, 0.0), (3.50E2, 0.025), (3.75E2, 0.1),
                      (3.94E2, 0.2), (4.00E2, 0.35)),
            'stress_columns': (0, ),
        },
        'section': {'type': 'HomogeneousSolidSection', 'thickness': None},
    },

    'Ti-6Al-4V': {
        'description': 'https://doi.org/10.1016/j.ijimpeng.2019.04.025',
        'density': 4.48E-9,
        'default_elastic': 'ISOTROPIC',
        'elastic': {
            'ISOTROPIC': {
                'type': 'ISOTROPIC',
                'table': ((1.287E5, 0.33), ),
                'stress_columns': (0, ),
            },
        },
        'plastic': {
            # (Yield stress, Plastic strain)
            'table': ((1.085E3, 0.00), (1.093E3, 0.02),
                      (1.123E3, 0.04), (1.155E3, 0.06)),
            'stress_columns': (0, ),
        },
        'section': {'type': 'HomogeneousSolidSection', 'thickness': None},
    },

    'Aluminum-7075': {
        'description': 'Aluminium alloy 7075, used for metallic ribs',
        'density': 2.81E-9,
        'default_elastic': 'ISOTROPIC',
        'elastic': {
            'ISOTROPIC': {
                'type': 'ISOTROPIC',
                'table': ((7.10E4, 0.33), ),
                'stress_columns': (0, ),
            },
        },
        'plastic': None,
        'section': {'type': 'HomogeneousSolidSection', 'thickness': None},
    },

    'IM7/8551-7': {
        'description': 'https://doi.org/10.1177/0021998312454478',
        'density': 1.272E-9,
        'default_elastic': 'LAMINA',
        'elastic': {
            'LAMINA': {
                'type': 'LAMINA',
                'table': ((165000.0, 8400.0, 0.34, 5600.0, 5600.0, 2800.0), ),
                'stress_columns': (0, 1, 3, 4, 5),
            },
            'ENGINEERING_CONSTANTS': {
                'type': 'ENGINEERING_CONSTANTS',
                'table': ((165000.0, 8400.0, 8400.0, 0.34, 0.34, 0.5,
                           5600.0, 5600.0, 2800.0), ),
                'stress_columns': (0, 1, 2, 6, 7, 8),
            },
        },
        'plastic': None,
        'section': {'type': 'HomogeneousSolidSection', 'thickness': None},

        #* Hashin damage, only applied when `pMesh['failure_model'] == 'hashin'`.
        #* Defined in the N-mm-tonne system only, see `hashin_damage_tables()`.
        #* *DAMAGE INITIATION, CRITERION=HASHIN: (XT, XC, YT, YC, SL, ST)
        #* *DAMAGE EVOLUTION, TYPE=ENERGY: (Gft, Gfc, Gmt, Gmc) in kJ/m^2 = N/mm
        'hashin': {
            'initiation_table': ((2560.0, 1590.0, 73.0, 185.0, 90.0, 92.5), ),
            'alpha': 1.0,
            'evolution_table': ((92.0, 80.0, 0.21, 0.8), ),
        },
    },
}


def get_material(name):
    '''
    Get the data of a material in `MATERIAL_LIBRARY`.

    Parameters
    --------------
    name: str
        name of the material, also the name used in the Abaqus model

    Returns
    --------------
    material: dict
        material data
    '''
    if name not in MATERIAL_LIBRARY:
        raise KeyError(
            'Unknown material "%s". Available: %s' % (name, list_materials()))

    return MATERIAL_LIBRARY[name]


def list_materials():
    '''
    Names of all materials in `MATERIAL_LIBRARY`.
    '''
    return sorted(MATERIAL_LIBRARY.keys())


def scale_table(table, stress_columns, unit_length='mm'):
    '''
    Convert a material table from the N-mm-tonne system to the target system.

    Parameters
    --------------
    table: tuple of tuples
        material table as stored in `MATERIAL_LIBRARY`

    stress_columns: tuple of int
        indices of the columns that carry a stress unit.
        The other columns are dimensionless (Poisson's ratio, plastic strain).

    unit_length: str
        'mm' or 'm', unit of length of the Abaqus model

    Returns
    --------------
    table: tuple of tuples
        converted material table
    '''
    if unit_length not in UNIT_SCALE:
        raise ValueError(
            'Wrong unit_length "%s", should be one of %s' % (
                unit_length, sorted(UNIT_SCALE.keys())))

    factor = UNIT_SCALE[unit_length]['stress']

    if factor == 1.0:
        return table

    return tuple(
        tuple(value * factor if i_col in stress_columns else value
              for i_col, value in enumerate(row))
        for row in table)


def scale_density(density, unit_length='mm'):
    '''
    Convert a density from tonne/mm^3 to the target unit system.
    '''
    if unit_length not in UNIT_SCALE:
        raise ValueError(
            'Wrong unit_length "%s", should be one of %s' % (
                unit_length, sorted(UNIT_SCALE.keys())))

    return density * UNIT_SCALE[unit_length]['density']
