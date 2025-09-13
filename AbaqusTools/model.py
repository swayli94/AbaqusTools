'''
Classes for Abaqus modeling procedures:

-   `Model`: class for creating an Abaqus model,
    including model initialization, building `Part` s, and procedures in
    *Assembly*, *Step*, *Interaction*, *Load* and *Job* modules;
    
-   `NodeOperation`: class for functions to manipulate nodes in an Abaqus model;

'''
import time
import copy

from AbaqusTools import IS_ABAQUS
from AbaqusTools.part import Part

if IS_ABAQUS:

    from abaqus import *
    from abaqusConstants import *
    from symbolicConstants import *
    from caeModules import *
    import mesh
    import odbAccess


class Model(object):
    '''
    Class for creating an Abaqus model:
    
    -   initialize the model, including attribute calculation
        and operations in the `Property` Abaqus Module;
    -   build parts with their `Part.build()` functions, including 
        operations in the `Sketch`, 'Part', and `Mesh` Abaqus Modules;
    -   carry out operations in the `Assembly`, `Step`, `Interaction`, 
        `Load` and `Job` Abaqus Modules; 
    
    Parameters
    --------------
    pGeo: dict, or None
        geometry parameters
        
    pMesh: dict, or None
        meshing parameters
        
    pRun: dict, or None
        load and running parameters
    
    Attributes
    --------------
    model: Abaqus model object
        Abaqus model
    '''
    def __init__(self, pGeo=None, pMesh=None, pRun=None, name_job='Job-1'):
        
        self.model = None
        
        self.pGeo = pGeo
        self.pMesh = pMesh
        self.pRun = pRun
        
        self.raise_exception = True
        
        self.name_job = str(name_job)

        self.initialization()
        
    @property
    def name_model(self):
        '''
        Abaqus default model name. DO not change it.
        '''
        return 'Model-1'

    def initialization(self):
        '''
        Initialize the Abaqus model:
        
        - create an Abaqus model
        - prepare attributes
        - carry out operations in the `Property` Abaqus Modules
        '''
        self.model = mdb.models[self.name_model]
        
        self.setup_property()

    #* =============================================
    #* Abaqus model property
    @property
    def rootAssembly(self):
        '''
        The `rootAssembly` of the model.
        '''
        return self.model.rootAssembly
    
    def instance(self, name_instance):
        '''
        The instance in the model with name `name_instance`.
        '''
        return self.model.rootAssembly.instances[name_instance]

    #* =============================================
    #* Abaqus procedure
    def build(self):
        '''
        Build an Abaqus model:
        
        - build parts with their `Part.build()` functions
        - operations in the `Assembly` Abaqus Modules
        - operations in the `Step` Abaqus Modules
        - operations in the `Interaction` Abaqus Modules
        - operations in the `Load` Abaqus Modules
        - operations in the `Job` Abaqus Modules
        
        '''
        t0 = time.time()
        
        self.setup_parts()
        self.setup_assembly()
        self.setup_steps()
        self.setup_outputs()
        self.setup_interactions()
        self.setup_loads()
        self.setup_jobs()

        t1 = time.time()

        print('>>> --------------------')
        print('    [Model] build time = %.1f (min)'%((t1-t0)/60.0))
        print('>>>')

    def setup_property(self):
        '''
        Define properties via Abaqus Module: Property
        
        - Create Material
        - Create Section
        '''

    def setup_parts(self):
        '''
        Build different parts with their `Part` objects.
        '''
    
    def setup_assembly(self):
        '''
        Assemble parts via Abaqus Module: Assembly
        
        - Create Instance
        - Linear Pattern
        - Radial Pattern
        - Translate Instance
        - Rotate Instance
        - Translate To
        - Create Constraint
        - Merge/Cut Instances
        '''
    
    def setup_steps(self):
        '''
        Define steps via Abaqus Module: Step
        
        - Create Step
        '''
    
    def setup_interactions(self):
        '''
        Define interactions via Abaqus Module: Interaction
        
        - Create Interaction
        - Create Interaction Property
        - Create Constraint
        - Create References Point
        '''
    
    def setup_loads(self):
        '''
        Define loads and boundary conditions via Abaqus Module: Load
        
        - Create Load
        - Create Boundary Condition
        - Create Predefined Field
        - Create Load Case
        '''
        
    def setup_outputs(self):
        '''
        Define outputs via Abaqus Module: Step
        
        - Create Field Output
        - Create History Output
        '''
        
    def setup_jobs(self):
        '''
        Define jobs via Abaqus Module: Job
        
        - Create Job
        '''
    
    #* =============================================
    #* Abaqus Assembly functions
    def create_mesh_node(self, name_part, x, y, z):
        '''
        Create a node in Assembly
        
        Parameters
        ----------------
        name_part: str
            name of the part that the nodes belong to.
            
        x, y, z: float
            coordinates of the node
        
        Return
        -----------------
        label: int
            node label
        '''
        node = self.model.parts[name_part].Node((x,y,z))
        label = node.label
        
        return label
    
    def create_geometry_set(self, name_set, name_instance, findAt_points, geometry='vertex', getClosest=True, searchTolerance=1000):
        '''
        Create a geometry set in Assembly
        
        Parameters
        ----------------
        name_set: str
            name of the set
            
        name_instance: str
            name of the instance that the geometry belong to.
            
        findAt_points: tuple, or a list of tuples
            point coordinates.
            
            Each face is specified by one point coordinate tuple (x,y,z).
            
            `findAt_points` can be either a tuple of one point, or a list of point tuples.
            
        geometry: str
            'vertex', 'edge', 'face', 'cell'
            
        getClosest: bool
            whether use `getClosest` to find an vertex or edge
            
        searchTolerance: float
            the distance within which the closest object must lie
        '''
        a = self.rootAssembly
        
        if geometry=='vertex':
            vertices = Part.get_vertices(self.instance(name_instance), findAt_points, getClosest, searchTolerance)
            a.Set(vertices=vertices, name=name_set)
        
        elif geometry=='edge':
            edges = Part.get_edges(self.instance(name_instance), findAt_points, getClosest, searchTolerance)
            a.Set(edges=edges, name=name_set)
        
        elif geometry=='face':
            faces = Part.get_faces(self.instance(name_instance), findAt_points)
            a.Set(faces=faces, name=name_set)
        
        elif geometry=='cell':
            cells = Part.get_cells(self.instance(name_instance), findAt_points)
            a.Set(cells=cells, name=name_set)
            
        else:
            raise Exception
            
    def create_node_set(self, name_set, name_instance, node_labels, unsorted=True):
        '''
        Create a node set in Assembly
        
        Parameters
        ----------------
        name_set: str
            name of the set
            
        name_instance: str
            name of the instance that the nodes belong to.
            
        node_labels: tuple of integers
            a tuple of node labels, e.g., (1,) or (1,2,3)
            
        unsorted: bool
            whether it is an unsorted node set
            
        References
        -----------------
        2.1.1 Node definition: Grouping nodes into node sets
        
            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/usb/default.htm?startat=pt01ch02s01aus05.html

            -   Node sets are used as convenient cross-references when defining loads, constraints, properties, etc. 

            -   By default, the nodes within a node set will be arranged in ascending order, and duplicate nodes will be removed. 
                Such a set is called a sorted node set. You may choose to create an unsorted node set as described later, 
                which is often useful for features that match two or more node sets. 
            
            -   For example, if you define multi-point constraints (General multi-point constraints, Section 28.2.2) between two node sets, 
                a constraint will be created between the first node in Set 1 and the first node in Set 2, 
                then between the second node in Set 1 and the second node in Set 2, etc. 
                It is important to ensure that the nodes are combined in the desired way. 
                Therefore, it is sometimes better to specify that a node set be stored in unsorted order.

        
        '''
        self.rootAssembly.SetFromNodeLabels(name=name_set, nodeLabels=((name_instance, node_labels),), unsorted=unsorted)

    def get_assembly_set(self, name_set, name_instance=None):
        '''
        Get a set from an instance
        
        Parameters
        -----------------
        name_set: str
            name of the set
            
        name_instance: str or None
        
            name of the instance that the set belongs to.
            
            If it is None, the set is directly defined in Assembly
            
        Return
        -----------------
        set: 
            Abaqus set object
            
        Examples
        ----------
        >>> set = mode.get_assembly_set(name_set, name_instance)
        '''
        if name_instance is None:
            return self.rootAssembly.sets['%s'%(name_set)]
        else:
            return self.rootAssembly.sets['%s.%s'%(name_instance, name_set)]

    def get_nodes_in_set(self, name_set, name_instance=None, return_labels=False):
        '''
        Get nodes in a set
        
        Parameters
        -----------------
        name_set: str
            name of the set of a face, edge, or vertex
            
        name_instance: str or None
        
            name of the instance that the set belongs to.
            
            If it is None, the set is directly defined in Assembly
            
        return_labels: bool
            whether also return labels of the nodes

        Return
        ---------------
        edge_nodes: MeshNodeArray object
            nodes on the edge
            
        label_nodes: list of integers
            label of nodes, optional output
            
        Examples
        ----------
        >>> nodes, label_nodes = model.get_nodes_in_set(name_set, name_instance, return_labels=True)
        >>> nodes = model.get_nodes_in_set(name_set, name_instance)
        '''
        set_ =self.get_assembly_set(name_set, name_instance)
        nodes = set_.nodes
        
        if return_labels:
            
            label_nodes = [node.label for node in nodes]
            return nodes, label_nodes
        
        else:
            
            return nodes

    #* =============================================
    #* Abaqus Load functions
    
    # The Analytical Field toolset
    # https://docs.software.vt.edu/abaqusv2022/English/SIMACAECAERefMap/simacae-m-Fld-sb.htm
    
    def create_analytical_mapped_field(self, name, xyzPointData, description='analytical mapped field',
                positiveNormalSearchTol=0.05, negativeNormalSearchTol=0.15, neighborhoodSearchTol=1000.0):
        '''
        Create an analytical mapped field
        
        Parameters
        -----------------
        name: str
            name of the analytical field
        
        description: str
            description of the analytical field
        
        xyzPointData: tuple of tuples
            a tuple of points, each point is a tuple of (x, y, z, value).
        
        positiveNormalSearchTol: float
            positive normal search distance tolerance
        
        negativeNormalSearchTol: float
            negative normal search distance tolerance
        
        neighborhoodSearchTol: float
            neighborhood search distance tolerance
        
        References
        -----------------
        Using analytical mapped fields
        https://docs.software.vt.edu/abaqusv2022/English/SIMACAECAERefMap/simacae-m-FldUsingmap-sb.htm
        
        About analytical mapped fields
        https://docs.software.vt.edu/abaqusv2022/English/SIMACAECAERefMap/simacae-c-fldusingmapapoint.htm#simacae-c-fldusingmapapoint
        
        Creating mapped fields from point cloud data
        https://docs.software.vt.edu/abaqusv2022/English/SIMACAECAERefMap/simacae-t-fldmapptcloudhelp.htm
        
        Point cloud data file formats for mapping.
        https://docs.software.vt.edu/abaqusv2022/English/?show=SIMACAECAERefMap/simacae-m-FldUsingmapPoint-sb.htm
        
        XYZ format
        https://docs.software.vt.edu/abaqusv2022/English/SIMACAECAERefMap/simacae-c-fldusingmappointxyz.htm
        '''
        self.model.MappedField(
            name=name, description=description, 
            regionType=POINT, partLevelData=False, localCsys=None, 
            pointDataFormat=XYZ, fieldDataType=SCALAR, 
            xyzPointData=xyzPointData, 
            positiveNormalSearchTol=positiveNormalSearchTol, 
            negativeNormalSearchTol=negativeNormalSearchTol, 
            neighborhoodSearchTol=neighborhoodSearchTol)

    def create_analytical_expression_field(self, name):
        '''
        Create an analytical expression field
        
        Parameters
        -----------------
        name: str
            name of the analytical field

        References
        -----------------
        Using analytical expression fields
        https://docs.software.vt.edu/abaqusv2022/English/SIMACAECAERefMap/simacae-m-FldUsing-sb.htm
        
        Building valid expressions
        https://docs.software.vt.edu/abaqusv2022/English/SIMACAECAERefMap/simacae-c-fldusingcreate.htm
        
        Creating expression fields
        https://docs.software.vt.edu/abaqusv2022/English/SIMACAECAERefMap/simacae-t-fldhelp.htm
        '''

    #* =============================================
    #* Abaqus Property functions
    def create_material_IM785517(self, elastic_type='LAMINA'):
        '''
        Create material: IM7/8551-7
        
        Parameters
        --------------
        elastic_type: str
            'LAMINA' or 'ENGINEERING_CONSTANTS'
        '''
        self.model.Material(name='IM7/8551-7', 
                            description='https://doi.org/10.1177/0021998312454478')
        
        if elastic_type == 'LAMINA':
            self.model.materials['IM7/8551-7'].Elastic(type=LAMINA, 
                table=((165000.0, 8400.0, 0.34, 5600.0, 5600.0, 2800.0), ))
        
        elif elastic_type == 'ENGINEERING_CONSTANTS':
            self.model.materials['IM7/8551-7'].Elastic(type=ENGINEERING_CONSTANTS, 
                table=((165000.0, 8400.0, 8400.0, 0.34, 0.34, 0.5, 5600.0, 5600.0, 2800.0), ))
        
        else:
            raise Exception
        
        
        if 'failure_model' not in self.pMesh:
            return
        
        elif self.pMesh['failure_model']=="LaRC05":
            '''
            When the model runs with LaRC05 user material (UMAT), user-defined output variables (UVARM).
            
            Reference:
             
                Miguel A. S. Matos , Silvestre T. Pinho, Failure criteria 
                for NCF composites Implementation in Abaqus,
                July 7, 2019, Imperial College London
                
                Chapter 7.2 UMAT
            '''   
            print('>>> --------------------')
            print('    [UMAT/UVARM] LaRC05')
            print('    Write the "PROPERTY TABLE" in the `*.inp` file.')
            print('>>>')
            
        elif self.pMesh['failure_model']=="Hashin":
            '''
            Abaqus tutorial:
            
                Use the following option to define the Hashin damage initiation criterion:
                *DAMAGE INITIATION, CRITERION=HASHIN, ALPHA, XT, XC, YT, YC, SL, ST
            
                https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/
                docs/v6.6/books/usb/default.htm?startat=pt05ch19s03abm41.html
                
                Use the following option to define the damage evolution law:
                *DAMAGE EVOLUTION, TYPE=ENERGY, SOFTENING=LINEAR, Gft, Gfc, Gmt, Gmc
                
                The 4 values in `table` are energies dissipated during damage for 
                fiber tension, fiber compression, matrix tension, 
                and matrix compression failure modes, respectively.

                https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/
                docs/v6.6/books/usb/default.htm?startat=pt05ch19s03abm42.html
            
            '''
            self.model.materials['IM7/8551-7'].Density(table=((1.272, ), ))
            
            self.model.materials['IM7/8551-7'].HashinDamageInitiation(table=((
                2560.0, 1590.0, 73.0, 185.0, 90.0, 92.5), ), alpha=1.0)
            
            self.model.materials['IM7/8551-7'].hashinDamageInitiation.DamageEvolution(
                type=ENERGY, table=((92.0, 80.0, 0.21, 0.8), ))

    def create_material_steel(self, unit_length='mm'):
        '''
        Create material: Steel
        
        Parameters
        --------------
        unit_length: str
            'm' or 'mm', unit of length in the Abaqus model
        '''
        self.model.Material(name='Steel')
        
        if unit_length == 'm':
        
            self.model.materials['Steel'].Density(table=((7.8E3, ), ))      # (Density (kg/m^3))
            self.model.materials['Steel'].Elastic(table=((2.1E11, 0.3),))   # (Young's modulus (N/m^2), Poisson ratio)
            self.model.materials['Steel'].Plastic(scaleStress=None, 
                table=((3.00E8, 0.0), (3.50E8, 0.025), (3.75E8, 0.1),       # (Yield stress (N/m^2), Plastic strain)
                       (3.94E8, 0.2), (4.00E8, 0.35)))
            
        elif unit_length == 'mm':
            
            self.model.materials['Steel'].Density(table=((7.8E-6, ), ))     # (Density (kg/mm^3))
            self.model.materials['Steel'].Elastic(table=((2.1E5, 0.3),))    # (Young's modulus (N/mm^2), Poisson ratio)
            self.model.materials['Steel'].Plastic(scaleStress=None, 
                table=((3.00E2, 0.0), (3.50E2, 0.025), (3.75E2, 0.1),       # (Yield stress (N/mm^2), Plastic strain)
                       (3.94E2, 0.2), (4.00E2, 0.35)))

        else:
            
            print('ERROR [create_material_steel]:')
            print('    Wrong unit_length input:', unit_length)
            print('    Should only be m or mm')
            
            raise Exception()

    def create_material_titanium(self, unit_length='mm'):
        '''
        Create material: Ti-6Al-4V
        
        https://doi.org/10.1016/j.ijimpeng.2019.04.025
        
        Parameters
        --------------
        unit_length: str
            'm' or 'mm', unit of length in the Abaqus model
        '''
        self.model.Material(name='Ti-6Al-4V')
        
        if unit_length == 'm':
        
            self.model.materials['Ti-6Al-4V'].Density(table=((4.48E3, ), ))     # (Density (kg/m^3))
            self.model.materials['Ti-6Al-4V'].Elastic(table=((1.287E11, 0.33),))# (Young's modulus (N/m^2), Poisson ratio)
            self.model.materials['Ti-6Al-4V'].Plastic(scaleStress=None, 
                table=((1.085E9, 0.00), (1.093E9, 0.02),                    # (Yield stress (N/m^2), Plastic strain)
                       (1.123E9, 0.04), (1.155E9, 0.06),))
            
        elif unit_length == 'mm':
            
            self.model.materials['Ti-6Al-4V'].Density(table=((4.48E-6, ), ))    # (Density (kg/mm^3))
            self.model.materials['Ti-6Al-4V'].Elastic(table=((1.287E5, 0.33),)) # (Young's modulus (N/mm^2), Poisson ratio)
            self.model.materials['Ti-6Al-4V'].Plastic(scaleStress=None, 
                table=((1.085E3, 0.00), (1.093E3, 0.02),                    # (Yield stress (N/mm^2), Plastic strain)
                       (1.123E3, 0.04), (1.155E3, 0.06),))

        else:
            
            print('ERROR [create_material_titanium]:')
            print('    Wrong unit_length input:', unit_length)
            print('    Should only be m or mm')
            
            raise Exception()

    def create_section_steel(self):
        '''
        Create section: Steel
        '''
        self.model.HomogeneousSolidSection(name='Steel', material='Steel', thickness=None)

    def create_section_titanium(self):
        '''
        Create section: Ti-6Al-4V
        '''
        self.model.HomogeneousSolidSection(name='Ti-6Al-4V', material='Ti-6Al-4V', thickness=None)

    def create_section_IM785517(self):
        '''
        Create section: IM7/8551-7
        '''
        self.model.HomogeneousSolidSection(name='IM7/8551-7', material='IM7/8551-7', thickness=None)

    def write_IM785517_property_table_inp(self, method='UMAT', fname_input='Job_1.inp'):
        '''
        Modify the `*.inp` file to add the "PROPERTY TABLE" for LaRC05.        
        Only works when the model runs with LaRC05's user material (UMAT).
        
        Parameters
        --------------
        method: str
            'UMAT' or 'UVARM'
        
        fname_input: str
            file name of the input
        
        References
        -----------------
        Miguel A. S. Matos , Silvestre T. Pinho, Failure criteria 
        for NCF composites Implementation in Abaqus,
        July 7, 2019, Imperial College London
        
        Chapter 7.2 UMAT
        '''
        if not self.pMesh['failure_model']=="LaRC05":
            return

        N_SKIP_LINE_IM78551 = 3
        
        if method == 'UMAT':
        
            PROPERTY_TABLE_IM78551 = [
                '*PROPERTY TABLE TYPE, NAME=LARC05PROPERTIES, PROPERTIES=19',
                '*MATERIAL, NAME=IM7/8551-7',
                '*TRANSVERSE SHEAR, TYPE=ORTHOTROPIC',
                '5.6E3, 5.6E3, 2.8E3',
                '*USER MATERIAL, CONSTANTS=0',
                '*DEPVAR',
                '9,',
                '1, FI_1_Matrix, Matrix failure index',
                '2, FI_2_Split, Matrix splitting failure index',
                '3, FI_3_FibreTen, Fibre tension failure index',
                '4, FI_4_Kink, Fibre kinking failure index',
                '5, FI_MAX, Maximum local failure index',
                '6, DMG_1, Damage for matrix failure',
                '7, DMG_2, Damage for matrix splitting',
                '8, DMG_3, Damage for tensile fibre failure',
                '9, DMG_4, Damage for fibre kinking',
                '*PROPERTY TABLE, TYPE=LARC05PROPERTIES',
                '1.6500E+5,  8.4000E+3,  3.4000E-1,  5.6000E+3,  2.5600E+3,  1.5900E+3,  7.3000E+1,  1.8500E+2, ',
                '9.0000E+1,  5.3000E+1,  8.2000E-2, -1.0000E+0,  6.0000E+1, -1.0000E+0,  2.1000E-1,  8.0000E-1, ',
                '9.2000E+1,  8.0000E+1,  1.2100E+0',
            ]
            
        elif method == 'UVARM':
            
            PROPERTY_TABLE_IM78551 = [
                '*PROPERTY TABLE TYPE, NAME=LARC05PROPERTIES, PROPERTIES=14',
                '*MATERIAL, NAME=IM7/8551-7',
                '*USER OUTPUT VARIABLES',
                '7,',
                '*ELASTIC, type=ENGINEERING CONSTANTS',
                '165000.0,8400.,8400., 0.34, 0.34, 0.5, 5600.0, 5600.0,',
                '2800.0,',
                '*PROPERTY TABLE, TYPE=LARC05PROPERTIES',
                '1.6500E+5,  8.4000E+3,  3.4000E-1,  5.6000E+3,  2.5600E+3,  1.5900E+3,  7.3000E+1,  1.8500E+2,',
                '9.0000E+1,  5.3000E+1,  8.2000E-2, -1.0000E+0, 90, 63',
            ]
            
        else:
            
            print('ERROR [write_IM785517_property_table_inp]:')
            print('      Unknown method input: ', method)
            print('      Only UMAT and UVARM are valid')
            raise Exception()

        with open(fname_input, 'r') as f:
            lines = f.readlines()

        with open(fname_input, 'w') as f:
            
            overwrite = False
            counter = 0
            
            for i_line in range(len(lines)):
                
                line = lines[i_line].split()
                
                if len(line)>=2:
                    if line[0]=='*Material,' and line[1]=='name=IM7/8551-7':
                        
                        print('>>> --------------------')
                        print('    [%s] Found MATERIALS (IM7/8551-7) in %s'%(method, fname_input))
                        print('    Overwrite the "PROPERTY TABLE"')
                        print('>>>')
                        
                        for i in range(len(PROPERTY_TABLE_IM78551)):
                            f.write(PROPERTY_TABLE_IM78551[i] + '\n')
                        
                        overwrite = True
                    
                if overwrite and counter < N_SKIP_LINE_IM78551:

                    counter += 1
                    continue
                    
                else:
                    
                    f.write(lines[i_line])

    #* =============================================
    #* Abaqus Step functions
    def create_static_step(self, timePeriod=1.0, maxNumInc=10000, 
                initialInc=0.01, minInc=1E-15, maxInc=0.1, nlgeom=False):
        '''
        Create a static analysis step
        '''
        if nlgeom:
            nlgeom = ON
        else:
            nlgeom = OFF

        self.model.StaticStep(
            name=           'Loading', 
            previous=       'Initial', 
            description=    'Static simulation',
            timePeriod=     timePeriod, 
            maxNumInc=      maxNumInc, 
            initialInc=     initialInc, 
            minInc=         minInc, 
            maxInc=         maxInc, 
            nlgeom=         nlgeom)
        
    def create_dynamic_step(self, nlgeom=False):
        '''
        Create a dynamic analysis step.
        
        When conducting dynamic simulations, the element type needs to be 
        selected from the 'explicit' elemLibrary.
        '''
        if nlgeom:
            nlgeom = ON
        else:
            nlgeom = OFF
        
        self.model.ExplicitDynamicsStep(
            name=               'Loading', 
            previous=           'Initial', 
            description=        'Dynamic (explicit) simulation', 
            nlgeom=             nlgeom, 
            improvedDtMethod=   ON)

    @staticmethod
    def write_static_step_inp(fname_input='Job_1.inp', 
                timePeriod=1.0, initialInc=0.01, minInc=1E-15, maxInc=0.1):
        '''
        Modify the `*.inp` file to overwrite the static step settings.
        
        Parameters
        --------------
        fname_input: str
            file name of the input
        '''
        N_SKIP_LINE_STATIC = 2
        
        with open(fname_input, 'r') as f:
            lines = f.readlines()

        with open(fname_input, 'w') as f:
            
            overwrite = False
            counter = 0
            
            for i_line in range(len(lines)):
                
                line = lines[i_line].split()
                
                if len(line)>=1:
                    if line[0]=='*Static':
                        
                        print('>>> --------------------')
                        print('    Overwrite the "*Static" in %s'%(fname_input))
                        print('>>>')
                        
                        f.write('*Static \n')
                        f.write('%.2E, %.2E, %.2E, %.2E \n'%(initialInc, timePeriod, minInc, maxInc))
                        
                        overwrite = True
                    
                if overwrite and counter < N_SKIP_LINE_STATIC:

                    counter += 1
                    continue
                    
                else:
                    
                    f.write(lines[i_line])

    @staticmethod
    def write_output_field_frequency_interval(fname_input='Job_1.inp', 
                frequency=1, numIntervals=0):
        '''
        Modify the `*.inp` file to overwrite the output field frequency settings.
        
        Parameters
        --------------
        fname_input: str
            file name of the input
        
        frequency: int
            output field in every n increments.
            Only output the last step when frequency <= 0.
        
        numIntervals: int
            output field in evenly spaced time intervals.
            The numIntervals sepcifies the number of outputs.
            Only take effects when it is larger than 0. 
        
        '''
        N_SKIP_LINE_STATIC = 1
        
        with open(fname_input, 'r') as f:
            lines = f.readlines()

        with open(fname_input, 'w') as f:
            
            overwrite = False
            counter = 0
            
            for i_line in range(len(lines)):
                
                line = lines[i_line].split()
                
                if len(line)>=2:
                    if line[0]=='*Output,' and (line[1]=='field' or line[1]=='field,'):
                        
                        print('>>> --------------------')
                        print('    Overwrite the "*Output, field" in %s'%(fname_input))
                        print('>>>')
                        
                        if numIntervals > 0:
                            f.write('*Output, field, number interval=%d \n'%(numIntervals))
                        
                        elif frequency <= 0:
                            f.write('*Output, field, frequency=99999 \n')
                            
                        elif frequency > 1:
                            f.write('*Output, field, frequency=%d \n'%(frequency))
                            
                        else:
                            f.write('*Output, field \n')

                        overwrite = True
                    
                if overwrite and counter < N_SKIP_LINE_STATIC:

                    counter += 1
                    continue
                    
                else:
                    
                    f.write(lines[i_line])

    #* =============================================
    #* Abaqus Interaction functions
    def create_interaction_property_contact(self):
        '''
        Create the interaction property of contact with friction
        '''
        name = 'IP-Contact-Friction'
        
        miu = 0.5
        if self.pRun is not None:
            if 'contact_friction_coef' in self.pRun.keys():
                miu = self.pRun['contact_friction_coef']
                if miu <= 0.0:
                    print('>>> [Contact property]:')
                    print('    No friction in tangential behaviors')

        # https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch22pyo19.html
        self.model.ContactProperty(name)

        # https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch22pyo19.html
        self.model.interactionProperties[name].TangentialBehavior(
            formulation=PENALTY,    # Specifies the friction coefficient formulation, e.g., PENALTY and EXPONENTIAL_DECAY.
            directionality=ISOTROPIC, slipRateDependency=OFF, 
            pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
            table=((miu, ), ),      # If formulation=PENALTY, the table data specify the following:
                                    # 1) friction coefficient in the slip direction, miu;
                                    # 2) slip rate, if the data depend on slip rate;
                                    # etc.
            shearStressLimit=None, maximumElasticSlip=FRACTION, 
            fraction=0.005, elasticSlipStiffness=None)
        
        # https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch22pyo19.html
        self.model.interactionProperties[name].NormalBehavior(
            pressureOverclosure=HARD, 
            allowSeparation=ON, 
            constraintEnforcementMethod=DEFAULT)

    def create_contact_constraints(self, name_face_pairs):
        '''
        Define general interaction for contacts, also define constraint `Tie` for face pairs. 
        These face pairs are excluded from the 'General Interaction'.
        
        Parameters
        ---------------
        name_face_pairs: list of tuples
        
            [(name of instance #1, name of main surface in #1, 
            name of instance #2, name of secondary surface in #2, name of the constraint), ...]
            
            The secondary surface should be smaller than the main surface,
            because the forbidden nodes for periodic boundary conditions should be 
            the nodes in the shorter edge of the secondary surface.
        '''
        self.model.ContactStd(name='General Interaction', createStepName='Initial')
        
        gi = self.model.interactions['General Interaction']

        gi.includedPairs.setValuesInStep(stepName='Initial', useAllstar=ON)

        for name_i1, name_s1, name_i2, name_s2, name_constraint in name_face_pairs:

            r1 = self.instance(name_i1).surfaces[name_s1]
            r2 = self.instance(name_i2).surfaces[name_s2]
        
            gi.excludedPairs.setValuesInStep(stepName='Initial', addPairs=((r1, r2), ))

            self.model.Tie(name='Tie-%s'%(name_constraint), main=r1, secondary=r2, 
                positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, 
                constraintEnforcement=NODE_TO_SURFACE, thickness=ON)

        gi.contactPropertyAssignments.appendInStep(stepName='Initial', 
                assignments=((GLOBAL, SELF, 'IP-Contact-Friction'), ))

    #* =============================================
    #* Abaqus Load functions
    def create_amplitude(self):
        '''
        Create amplitude functions
        '''
        self.model.TabularAmplitude(name='Constant-Amp', timeSpan=STEP,
                smooth=SOLVER_DEFAULT, data=((0.0, 1.0), (1.0, 1.0)))
        
        self.model.TabularAmplitude(name='Ramp-Amp', timeSpan=STEP,
                smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
    
    def create_reference_point(self, x, y, z, name_rp):
        '''
        Create a reference point feature in Assembly by coordinates, and specify its name
        
        Parameters
        ------------
        name_rp: str
            name of the reference point
        '''
        self.rootAssembly.ReferencePoint(point=(x, y, z))
        self.rootAssembly.features.changeKey(fromName='RP-1', toName=name_rp)
    
    def get_reference_point(self, name_rp):
        '''
        Get a reference point feature by name
        
        Parameters
        ---------------
        name_rp: str
            name of the reference point
        
        Return
        ---------------
        rp: ReferencePoint object
            
        Examples
        ----------
        >>> rp = model.get_reference_point(name_rp)
        '''
        id_rp = self.rootAssembly.features[name_rp].id
        return self.rootAssembly.referencePoints[id_rp]
    
    def create_reference_point_set(self, name_set, name_rp):
        '''
        Create a geometry set for a reference point
        
        Parameters
        ------------
        name_set: str
            name of the set
            
        name_rp: str
            name of the reference point
        '''
        rp = self.get_reference_point(name_rp)
        self.rootAssembly.Set(referencePoints=(rp,), name=name_set)
    
    #* =============================================
    #* Abaqus Job functions
    def submit_job(self, name_job=None, only_data_check=False):
        '''
        Submit job
        
        Parameters
        ---------------
        name_job: str, None
            name of the job.
            Default is None, which uses the default job name
            
        only_data_check: bool
            whether only carries out data check
        '''
        t0 = time.time()
        
        if name_job is None:
            name_job = self.name_job
        
        mdb.jobs[name_job].submit(consistencyChecking=OFF, datacheckJob=only_data_check)

        t1 = time.time()

        print('>>> --------------------')
        print('    [Job] %s, time = %.1f (min)'%(name_job, (t1-t0)/60.0))
        print('>>>')

    def write_job_inp(self, name_job=None):
        '''
        Write the job input file
        
        Parameters
        ---------------
        name_job: str, None
            name of the job.
            Default is None, which uses the default job name
        '''
        if name_job is None:
            name_job = self.name_job
        
        mdb.jobs[name_job].writeInput(consistencyChecking=OFF)

    #* =============================================
    #* Abaqus other functions
    def set_view(self):
        '''
        Set Abaqus view
        '''
        a = self.rootAssembly
        a.regenerate()
        viewport = session.viewports['Viewport: 1']
        viewport.setValues(displayedObject=a)
        viewport.assemblyDisplay.setValues(optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
        viewport.assemblyDisplay.geometryOptions.setValues(datumPoints=OFF, datumAxes=OFF, datumPlanes=OFF, datumCoordSystems=OFF)
        viewport.setColor(globalTranslucency=True)
        viewport.view.rotate(xAngle=20, yAngle=60, zAngle=10, mode=TOTAL)
        session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='#FFFFFF')
        session.viewports["Viewport: 1"].viewportAnnotationOptions.setValues(compass=OFF)

    def set_view_fixed_origin(self):
        '''
        Set Abaqus view where the location of the origin (0,0,0) is fixed.
        '''
        session.View(name='User-1', 
            nearPlane=2549.1, farPlane=4259.2, 
            width=1669.3, height=809, 
            projection=PERSPECTIVE, 
            cameraPosition=(-1797.7, 1551.9, 2152.4), 
            cameraUpVector=(0.42948, 0.88149, -0.19628), 
            cameraTarget=(771.82, -44.15, 607.24), 
            viewOffsetX=-3.0688, viewOffsetY=-3.6258, autoFit=OFF)
        
        session.View(name='User-2', 
            nearPlane=2962.9, farPlane=5321.3, 
            width=2342.8, height=1207.5, 
            projection=PERSPECTIVE, 
            cameraPosition=(-2122.3, 2246.5, 3284.7), 
            cameraUpVector=(0.51875, 0.82747, -0.21495), 
            cameraTarget=(839.78, -81.679, 1470.9), 
            viewOffsetX=38.909, viewOffsetY=-64.822, autoFit=OFF)
        
        session.viewports['Viewport: 1'].view.setValues(session.views['User-1'])

    @staticmethod
    def save_cae(pathName):
        '''
        Saves an Mdb object to disk at the specified location
        '''
        mdb.saveAs(pathName=pathName)


class NodeOperation(object):
    '''
    A class that contains functions to manipulate nodes in Abaqus
    '''
    def __init__(self):
        pass

    @staticmethod
    def get_set(myMdl, name_set, name_part=None, name_instance=None):
        '''
        Get nodes on a face set.
        
        Parameters
        -----------------
        myMdl: Abaqus Model object
            model object
            
        name_set: str
        
            name of the face set.

        name_part: str or None
        
            name of the part that the set belongs to.
            
            If it is None, the set is defined in Assembly
            
        name_instance: str or None
        
            name of the instance that the set belongs to.
            
            If it is None, the set is directly defined in Assembly
        
        Return
        ---------------
        set: Set object
            Abaqus set
        '''
        if name_part is not None:
            set_= myMdl.parts[name_part].sets['%s'%(name_set)]
        elif name_instance is None:
            set_= myMdl.rootAssembly.sets['%s'%(name_set)]
        else:
            set_= myMdl.rootAssembly.sets['%s.%s'%(name_instance, name_set)]
        
        return set_
        
    @staticmethod
    def get_nodes_from_face(myMdl, name_face_set, coordinates_for_sorting, name_part=None, name_instance=None):
        '''
        Get nodes on a face set.
        
        Parameters
        -----------------
        myMdl: Abaqus Model object
            model object
            
        name_face_set: str
        
            name of the face set.

        coordinates_for_sorting: tuple of int
        
            a tuple contains coordinate indices for sorting, e.g., (0,1), (1,2), (2,0).
            
            When the face pairs have the same node distribution,
            the node pairs can be found by sorting the nodes in each face in the same manner.
        
        name_part: str or None
        
            name of the part that the set belongs to.
            
            If it is None, the set is defined in Assembly
            
        name_instance: str or None
        
            name of the instance that the set belongs to.
            
            If it is None, the set is directly defined in Assembly
        
        Return
        ---------------
        face_nodes: MeshNodeArray object
            nodes on the face
        '''
        set_ = NodeOperation.get_set(myMdl, name_face_set, name_part, name_instance)
        
        #* Sort nodes by coordinates to find pairs automatically
        if len(coordinates_for_sorting)==2:

            coordA, coordB = coordinates_for_sorting
            face_nodes = sorted(set_.nodes, key=lambda node: (node.coordinates[coordA], node.coordinates[coordB]))
            
        elif len(coordinates_for_sorting)==1:
        
            coord = coordinates_for_sorting[0]
            face_nodes = sorted(set_.nodes, key=lambda node: node.coordinates[coord])
        
        else:
            
            print('>>> --------------------')
            print('[Error]: NodeOperation.get_nodes_from_face')
            print('         The `coordinates_for_sorting` should be a tuple of one integer for 2D models or two integers for 3D model,')
            print('         but it is ', coordinates_for_sorting)
            raise Exception()

        return face_nodes

    @staticmethod
    def exclude_forbidden_nodes(myMdl, nodes, name_forbidden_sets=[], label_forbidden_nodes=[], name_part=None, name_instance=None):
        '''
        Exclude forbidden nodes from an array of nodes, and exclude the corresponding node from the slave face.
        
        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        nodes: MeshNodeArray object
            nodes to be selected
            
        name_forbidden_sets: list of strings
        
            name of sets that consists of forbidden nodes.
            
            It can be a node set, or a geometry set of vertices, edges, or faces.
            
        label_forbidden_nodes: list of integers
        
            labels of forbidden nodes.
            
            An overlap between `name_forbidden_sets` and `label_forbidden_nodes` is allowed.
            
        name_part: str or None
        
            name of the part that the set belongs to.
            
            If it is None, the set is defined in Assembly
            
        name_instance: str or None
        
            name of the instance that the set belongs to.
            
            If it is None, the set is directly defined in Assembly
            
        Return
        ---------------
        label_nodes: list of integers
            labels of the remaining nodes
            
        label_forbidden: list of integers
            labels of forbidden nodes
        
        Examples
        ---------------
        >>> label_nodes, label_forbidden = \\
        >>>         NodeOperation.exclude_forbidden_nodes(myMdl, nodes, 
        >>>         name_forbidden_sets=[], label_forbidden_nodes=[],
        >>>         name_part=None, name_instance=None)
        '''
        label_forbidden = copy.deepcopy(label_forbidden_nodes)
        label_nodes = []

        for name_set in name_forbidden_sets:
            
            set_ = NodeOperation.get_set(myMdl, name_set, name_part, name_instance)
            label_forbidden += [node.label for node in set_.nodes]

        for node in nodes:

            # Check master node is not among the forbidden nodes list
            if node.label in label_forbidden:
                continue
            
            label_nodes.append(node.label)
            
            # Update forbidden nodes list
            label_forbidden.append(node.label)
                
        return label_nodes, label_forbidden

    @staticmethod
    def create_node_set(myMdl, nodes, name_node_set, name_part=None, name_instance=None):
        '''
        Create an unsorted node set
        
        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        nodes: MeshNodeArray object
            nodes, e.g., sets['*'].nodes

        name_set: str
            name of the node set
            
        name_part: str or None
        
            name of the part that the set belongs to.
            
            One of name_part and name_instance must be provided.
            
        name_instance: str or None
        
            name of the instance that the set belongs to.
            
            One of name_part and name_instance must be provided.
        '''
        label_nodes = [node.label for node in nodes]
        
        if name_part is not None and name_instance is None:

            myMdl.parts[name_part].SetFromNodeLabels(name=name_node_set, nodeLabels=label_nodes, unsorted=True)
            
        elif name_instance is not None and name_part is None:
            myMdl.rootAssembly.SetFromNodeLabels(name=name_node_set, nodeLabels=((name_instance, label_nodes),), unsorted=True)
            
        else:
            print('>>> Error [NodeOperation.create_node_set]')
            print('    Only one of name_part and name_instance must be provided.')
            print('name_part     ', name_part)
            print('name_instance ', name_instance)
            raise Exception()
        
    @staticmethod
    def create_face_node_set(myMdl, name_face_set, coords_sorting, name_forbidden_sets=[], label_forbidden_nodes=[], 
                                name_part=None, name_instance=None, name_node_set=None):
        '''
        Create node sets on the master/slave faces for the periodic boundary condition.
        
        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        coords_sorting: tuple of int
            a tuple contains coordinate indices for sorting nodes, e.g., (0,1), (1,2), (2,0).
            
        name_forbidden_sets: list of strings
        
            name of sets that consists of forbidden nodes.
            
            It can be a node set, or a geometry set of vertices, edges, or faces.
            
        label_forbidden_nodes: list of integers
        
            labels of forbidden nodes.
            
            An overlap between `name_forbidden_sets` and `label_forbidden_nodes` is allowed.
            
        name_part: str or None
        
            name of the part that the set belongs to.
            
            One of name_part and name_instance must be provided.
            
        name_instance: str or None
        
            name of the instance that the set belongs to.
            
            One of name_part and name_instance must be provided.
            
        name_node_set: str, or None
            name of the node set. If None, use default name.
            
        Returns
        -----------------
        name_node_set: str
            name of the node set
        '''
        #* Sort nodes by coordinates to find node pairs automatically
        face_nodes = NodeOperation.get_nodes_from_face(myMdl, name_face_set, coords_sorting, name_part, name_instance)

        #* Exclude forbidden nodes
        label_nodes, label_forbidden_nodes = \
            NodeOperation.exclude_forbidden_nodes(myMdl, face_nodes, name_forbidden_sets, label_forbidden_nodes,
                                                    name_part, name_instance)
        if name_node_set is None:
            name_node_set = 'FNode-%s-%s'%(name_instance, name_face_set)

        #* Create node sets
        if name_part is not None and name_instance is None:

            myMdl.parts[name_part].SetFromNodeLabels(name=name_node_set, nodeLabels=label_nodes, unsorted=True)
            
        elif name_instance is not None and name_part is None:
            myMdl.rootAssembly.SetFromNodeLabels(name=name_node_set, nodeLabels=((name_instance, label_nodes),), unsorted=True)
            
        else:
            print('>>> Error [NodeOperation.create_face_node_set]')
            print('    Only one of name_part and name_instance must be provided.')
            print('name_part     ', name_part)
            print('name_instance ', name_instance)
            raise Exception()
        
        return name_node_set

