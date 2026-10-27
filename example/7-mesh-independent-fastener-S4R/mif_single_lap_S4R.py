'''
Single lap specimen for mesh-independent fasteners (S4R).
'''
import os
import time
import numpy as np
import json
from AbaqusTools import Part, IS_ABAQUS

if IS_ABAQUS:
    from abaqus import *
    from abaqusConstants import *
    from sketch import *
    import mesh
    
from AbaqusTools import Part, Model
from AbaqusTools import OdbOperation


class ImplicitHolePlate(Part):
    '''
    Plate with an implicit hole.
    '''
    def __init__(self, model, name_part, pMesh,
            length_x, length_y, length_z, xc_hole, yc_hole, r_hole):
        
        super(ImplicitHolePlate, self).__init__(model, pGeo=None, pMesh=pMesh)
        
        self.name_part = name_part

        #* Attributes
        self.length_x = length_x
        self.length_y = length_y
        self.length_z = length_z

        self.xc_hole = xc_hole
        self.yc_hole = yc_hole
        self.r_hole = r_hole
        
        self.width_partition =self.pMesh['ImplicitModelling']['width_partition']
                
    def build(self):
        '''
        Build an Abaqus part:
        
        - operations in the `Sketch` Abaqus Modules
        
        - operations in the `Part` Abaqus Modules
        
            - Create Part
            - Create Partition
            - Create surfaces and sets
            
        - operations in the `Mesh` Abaqus Modules
        
            - Set seeding
            - Create Mesh
            - Assign Element Type
            
        - operations in the `Property` Abaqus Modules
        
            - Assign Section
            - Create Composite Layup
        '''
        t0 = time.time()
        
        #* Abaqus Module: Sketch
        self.create_sketch()
        
        #* Abaqus Module: Part
        self.create_part()
        self.create_partition()
        self.create_surface()
        self.create_set()
        
        self.create_partition_hole()

        if not self.is_only_geometry:

            #* Abaqus Module: Mesh
            self.set_seeding()
            self.create_mesh()
            self.set_element_type()

            #* Abaqus Module: Property
            self.set_section_assignment()
            self.set_composite_layups()
        
        t1 = time.time()

        print('>>> --------------------')
        print('    [Part: %s] build time = %.1f (min)'%(self.name_part, (t1-t0)/60.0))
        print('>>>')

    def create_sketch(self):
        '''
        Create the sketch in X-Y plane.
        '''
        #* Sketch points in X-Y plane (a rectangle)
        sketch_points = np.zeros((4,2))
        sketch_points[0,:] = [0.0,        0.0]
        sketch_points[1,:] = [self.length_x, 0.0]
        sketch_points[2,:] = [self.length_x, self.length_y]
        sketch_points[3,:] = [0.0,        self.length_y]

        #* Create the open-hole plate sketch in x-y plane
        mySkt = self.model.ConstrainedSketch(name='plate_top_view', sheetSize=200)
        
        for i in range(4):
            mySkt.Line(point1=tuple(sketch_points[i-1,:]), 
                       point2=tuple(sketch_points[i,:]))

    def create_part(self):
        '''
        Create the part by extruding the sketch (Y-Z plane) in X direction.
        '''
        #* Create part
        myPrt = self.model.Part(name=self.name_part, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    
        #* Reference plane and axis
        myPrt.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
        self.rename_feature(myPrt, 'XYPLANE')
        
        myPrt.DatumAxisByPrincipalAxis(principalAxis=XAXIS)
        self.rename_feature(myPrt, 'XAXIS')
        
        myPrt.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
        self.rename_feature(myPrt, 'YAXIS')
        
        myPrt.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
        self.rename_feature(myPrt, 'ZAXIS')
        
        self.create_datum_csys_3p(myPrt, 'csys_plate', origin=[0.0, 0.0, 0.0],
                                    dx=[1, 0, 0], dy=[0, 1, 0])
        
        #* Plane for plate sketch
        transform = myPrt.MakeSketchTransform(
            sketchPlane=self.get_datum_by_name(myPrt, 'XYPLANE'),
            sketchUpEdge=self.get_datum_by_name(myPrt, 'XAXIS'), 
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 0.0))
    
        #* Section sketch
        mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=200, transform=transform)
        mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)
        mySkt.retrieveSketch(sketch=self.model.sketches['plate_top_view'])

        #* Part by sketch
        myPrt.BaseShell(sketch=mySkt)

        #* Post procedure
        myPrt.setValues(geometryRefinement=EXTRA_FINE)
        del self.model.sketches['__profile__']
    
    #* Surface, set for the entire part
    
    def create_partition(self):
        pass
     
    def create_surface(self):
        myPrt = self.model.parts[self.name_part]
        faces = self.get_faces(myPrt, (0.5*self.length_x, 0.5*self.length_y, 0.0),
                getClosest=True, searchTolerance=1E-3)
        myPrt.Surface(side1Faces=faces, name='all')

    def create_set(self):

        lx = self.length_x
        ly = self.length_y

        myPrt = self.model.parts[self.name_part]
        myPrt.Set(faces=myPrt.faces, name='all')
        
        self.create_geometry_set('edge_x0', (0.0,    0.5*ly, 0.0), geometry='edge')
        self.create_geometry_set('edge_x1', (lx,     0.5*ly, 0.0), geometry='edge')
        self.create_geometry_set('edge_y0', (0.5*lx, 0.0,    0.0), geometry='edge')
        self.create_geometry_set('edge_y1', (0.5*lx, ly,     0.0), geometry='edge')

        self.create_geometry_set('vertex_00', (0.0, 0.0, 0.0), geometry='vertex')
        self.create_geometry_set('vertex_10', (lx,  0.0, 0.0), geometry='vertex')
        self.create_geometry_set('vertex_01', (0.0, ly,  0.0), geometry='vertex')
        self.create_geometry_set('vertex_11', (lx,  ly,  0.0), geometry='vertex')

    #* Partition and create surfaces and sets for the partition
    
    def create_partition_hole(self):    
        '''
        After `create_surface` and `create_set`,
        partition a circle and square for the structure mesh around hole.
        '''        
        myPrt = self.model.parts[self.name_part]
        
        #* Partition cell to squares by 4 planes
        x0 = self.xc_hole - 0.5*self.width_partition
        x1 = self.xc_hole + 0.5*self.width_partition
        y0 = self.yc_hole - 0.5*self.width_partition
        y1 = self.yc_hole + 0.5*self.width_partition
        
        myPrt.PartitionFaceByShortestPath(faces=myPrt.faces, 
            point1=(x0, 0.0, 0.0), point2=(x0, self.length_y, 0.0))
        myPrt.PartitionFaceByShortestPath(faces=myPrt.faces, 
            point1=(x1, 0.0, 0.0), point2=(x1, self.length_y, 0.0))
        myPrt.PartitionFaceByShortestPath(faces=myPrt.faces, 
            point1=(0.0, y0, 0.0), point2=(self.length_x, y0, 0.0))
        myPrt.PartitionFaceByShortestPath(faces=myPrt.faces, 
            point1=(0.0, y1, 0.0), point2=(self.length_x, y1, 0.0))

        dd = 0.5*(self.width_partition*0.5+self.r_hole)
        pt = (self.xc_hole+dd, self.yc_hole, 0.0)
        self.create_geometry_set('partition_square', pt, geometry='face')
        
        faces = self.get_faces(myPrt, pt, getClosest=True, searchTolerance=1E-3)
        myPrt.Surface(side1Faces=faces, name='partition_square')
        
        self.create_geometry_set('edge_partition_x0', (x0, 0.5*(y0+y1), 0.0), geometry='edge')
        self.create_geometry_set('edge_partition_x1', (x1, 0.5*(y0+y1), 0.0), geometry='edge')
        self.create_geometry_set('edge_partition_y0', (0.5*(x0+x1), y0, 0.0), geometry='edge')
        self.create_geometry_set('edge_partition_y1', (0.5*(x0+x1), y1, 0.0), geometry='edge')
        
        myPrt.SetByBoolean(
            name='remainder',
            sets=(myPrt.sets['all'], myPrt.sets['partition_square']),
            operation=DIFFERENCE,
        )
        
    #* Meshing
    
    def set_seeding(self):

        myPrt = self.model.parts[self.name_part]
        myPrt.seedPart(size=self.pMesh['plate_seedPart_size'], 
                        deviationFactor=0.1, minSizeFactor=0.1)
        
        num = self.pMesh['fastener_partition_seedEdgeByNumber']
        
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_partition_x0'].edges, number=num, constraint=FIXED)
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_partition_x1'].edges, number=num, constraint=FIXED)
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_partition_y0'].edges, number=num, constraint=FIXED)
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_partition_y1'].edges, number=num, constraint=FIXED)

    def create_mesh(self):
        
        myPrt = self.model.parts[self.name_part]
        myPrt.setMeshControls(regions=myPrt.faces, elemShape=QUAD, technique=STRUCTURED)
        myPrt.generateMesh()

    def set_element_type(self):
        
        myPrt = self.model.parts[self.name_part]
        self.set_element_type_of_part(myPrt, kind='shell')
        
    def set_section_assignment(self):
        
        myPrt = self.model.parts[self.name_part]
        
        self.set_CompositeLayup_of_set(myPrt, 
                name_set=       'remainder', 
                total_thickness=self.length_z, 
                ply_angle=      self.pMesh['plate_CompositePly_orientationValue'],
                eNum_thickness= self.pMesh['num_element_thickness'],
                symmetric=      self.pMesh['plate_CompositeLayup_symmetric'],
                numIntPoints=   self.pMesh['plate_CompositePly_numIntPoints'],
                name_csys_datum='csys_plate',
                material='orthotropic',
                elementType='shell')

        self.set_CompositeLayup_of_set(myPrt, 
                name_set=       'partition_square', 
                total_thickness=self.length_z, 
                ply_angle=      self.pMesh['plate_CompositePly_orientationValue'],
                eNum_thickness= self.pMesh['num_element_thickness'],
                symmetric=      self.pMesh['plate_CompositeLayup_symmetric'],
                numIntPoints=   self.pMesh['plate_CompositePly_numIntPoints'],
                name_csys_datum='csys_plate',
                material='implicit_modelling',
                elementType='shell')


class SingleLapBoltedJoint(Model):
    
    def __init__(self, name_job, pGeo, pMesh, pRun, displacement=[1.0, 0.0, 0.0]):
        
        super(SingleLapBoltedJoint,self).__init__(pGeo=pGeo, pMesh=pMesh, pRun=pRun, name_job=name_job)
        
        self.displacement = displacement
        self.label_rp = 'RF_load'
        
    def initialization(self):
        
        self.model = mdb.models[str(self.name_model)]

        self.model.Material(name='orthotropic', 
                            description='custom orthotropic material with engineering constants')
        
        # Material properties (N-mm-MPa)
        e11 = self.pMesh['E11']
        e22 = self.pMesh['E22']
        e33 = self.pMesh['E33']
        nu12 = self.pMesh['nu12']
        nu13 = self.pMesh['nu13']
        nu23 = self.pMesh['nu23']
        g12 = self.pMesh['G12']
        g13 = self.pMesh['G13']
        g23 = self.pMesh['G23']

        self.model.materials['orthotropic'].Elastic(type=ENGINEERING_CONSTANTS, 
                table=((e11, e22, e33, nu12, nu13, nu23, g12, g13, g23), ))
        
        self.model.HomogeneousSolidSection(name='orthotropic', material='orthotropic', thickness=None)
        
        self._create_implicit_modelling_material()
        
    def _create_implicit_modelling_material(self):
        '''
        Create a material with implicit modelling properties.
        '''
        self.model.Material(name='implicit_modelling', 
                            description='material for implicit modelling')
        
        mat = self.pMesh['ImplicitModelling']
        
        # Material properties (N-mm-MPa)
        e11 = mat['E11']
        e22 = mat['E22']
        e33 = mat['E33']
        nu12 = mat['nu12']
        nu13 = mat['nu13']
        nu23 = mat['nu23']
        g12 = mat['G12']
        g13 = mat['G13']
        g23 = mat['G23']
        
        self.model.materials['implicit_modelling'].Elastic(type=ENGINEERING_CONSTANTS, 
                table=((e11, e22, e33, nu12, nu13, nu23, g12, g13, g23), ))
        
        self.model.HomogeneousSolidSection(name='implicit_modelling',
                                material='implicit_modelling', thickness=None)

    def setup_parts(self):
        
        self.length_x = self.pGeo['len_x_plate']
        self.length_y = self.pGeo['len_y_plate']
        self.length_z = self.pGeo['len_z_plate']
        self.r_hole = self.pGeo['fasteners'][0]['r_hole']
        self.r_washer = self.pGeo['fasteners'][0]['r_washer']
        self.thick_head = self.pGeo['fasteners'][0]['thick_head']
        self.xc_hole = self.pGeo['fasteners'][0]['x_center']
        self.yc_hole = self.pGeo['fasteners'][0]['y_center']
        
        self.plate_0 = ImplicitHolePlate(self.model,
                            name_part='plate_0',
                            pMesh=self.pMesh,
                            length_x=self.length_x,
                            length_y=self.length_y,
                            length_z=self.length_z,
                            xc_hole=self.xc_hole,
                            yc_hole=self.yc_hole,
                            r_hole=self.r_hole)
        self.plate_0.build()
        
        self.plate_1 = ImplicitHolePlate(self.model,
                            name_part='plate_1',
                            pMesh=self.pMesh,
                            length_x=self.length_x,
                            length_y=self.length_y,
                            length_z=self.length_z,
                            xc_hole=self.length_x-self.xc_hole,
                            yc_hole=self.yc_hole,
                            r_hole=self.r_hole)
        self.plate_1.build()
        
    def setup_assembly(self):
        a = self.rootAssembly
        a.Instance(name='plate_0', part=self.model.parts['plate_0'], dependent=ON)
        a.Instance(name='plate_1', part=self.model.parts['plate_1'], dependent=ON)

        dz_plate0 = 0.5*self.plate_0.length_z
        a.translate(instanceList=('plate_0',), vector=(0.0, 0.0, dz_plate0))

        dx_plate1 = self.xc_hole - (self.length_x - self.xc_hole)
        dz_plate1 = 0.5*self.plate_1.length_z + self.plate_0.length_z
        a.translate(instanceList=('plate_1',),
                    vector=(dx_plate1, 0.0, dz_plate1))

        hole_center_0 = [
            self.plate_0.xc_hole,
            self.plate_0.yc_hole,
            dz_plate0
        ]

        hole_center_1 = [
            self.plate_1.xc_hole + dx_plate1,
            self.plate_1.yc_hole,
            dz_plate1
        ]

        self.create_reference_point(hole_center_0[0], hole_center_0[1], hole_center_0[2], 'RP_hole_0')
        self.create_reference_point_set('RP_hole_0', 'RP_hole_0')  

        self.create_reference_point(hole_center_1[0], hole_center_1[1], hole_center_1[2], 'RP_hole_1')
        self.create_reference_point_set('RP_hole_1', 'RP_hole_1')  

    def setup_steps(self):
        self.create_static_step(
            timePeriod= self.pRun['timePeriod'],
            maxNumInc=  self.pRun['maxNumInc'],
            initialInc= self.pRun['initialInc'],
            minInc=     self.pRun['minInc'],
            maxInc=     self.pRun['maxInc'],
            nlgeom=     self.pRun['nlgeom'],
        )
        
    def setup_interactions(self):

        self._setup_mesh_independent_fastener_interactions()
        
    def _setup_mesh_independent_fastener_interactions(self):
        '''
        Define interactions for mesh-independent fastener.
        '''
        a = self.rootAssembly
        r_washer = self.pGeo['fasteners'][0]['r_washer']
        r_hole = self.pGeo['fasteners'][0]['r_hole']

        self.model.ConnectorSection(name='ConnSect-0', assembledType=BEAM)

        targetSurfaces = (
            a.instances['plate_0'].surfaces['partition_square'],
            a.instances['plate_1'].surfaces['partition_square'],
        )

        a.engineeringFeatures.PointFastener(name='Fastener-0',
            region=a.sets['RP_hole_0'],
            targetSurfaces=targetSurfaces, 
            directionVector=((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
            maximumLayers=len(targetSurfaces), 
            weightingMethod=LINEAR,
            physicalRadius=r_hole,
            sectionName='ConnSect-0',
            connectionType=CONNECTOR, 
            unsorted=OFF)
        
    def setup_loads(self):

        a = self.rootAssembly
        
        self.create_reference_point(self.length_x*2, 0.5*self.length_y, 0, self.label_rp)
        self.create_reference_point_set(self.label_rp, self.label_rp)  

        #* Encastre BC on the left edge of plate_0
        self.model.EncastreBC(name='BC-x0', createStepName='Initial',
            region=a.sets['plate_0.edge_x0'],
            localCsys=None)
        
        #* Coupling BC on the right edge of plate_1
        # Only constrain the loaded DOFs so the end face can deform freely in other directions
        _d = self.displacement
        self.model.Coupling(name='Coupling-plate_1_edge_x1',
            controlPoint=a.sets[self.label_rp],
            surface=a.instances['plate_1'].sets['edge_x1'],
            influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
            alpha=0.0, localCsys=None,
            u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
        
        # #* Apply displacement BCs on the reference points
        self.model.DisplacementBC(name=self.label_rp, createStepName='Loading', 
            region=a.sets[self.label_rp],
            u1=_d[0], u2=_d[1], u3=_d[2],
            ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF,
            distributionType=UNIFORM, fieldName='', localCsys=None)
        
    def setup_outputs(self):

        self.model.FieldOutputRequest(name='F-Output-1', 
            createStepName='Loading', variables=('S', 'E', 'SE', 'SF', 'U', 'RF'),
            frequency=LAST_INCREMENT)
        
        # SF: section forces and moments (N11, N22, N12, M11, M22, M12)
        # N**: in-plane forces per unit length (N/mm)
        # M**: bending moments per unit length (N*mm/mm)
        # SE: section strains of mid-plane (epsilon11, epsilon22, epsilon12, kappa11, kappa22, kappa12)
        
        variables = ('S', 'TSHR', 'E')
        if 'failure_model' in self.pMesh:
            if self.pMesh['failure_model'] == 'Hashin':
                variables_hashin = ('DMICRT', 'HSNFTCRT', 'HSNFCCRT', 'HSNMTCRT', 'HSNMCCRT')
                variables += variables_hashin

        self.model.FieldOutputRequest(name='Layup-Output-0', 
            createStepName='Loading', variables=variables,
            frequency=LAST_INCREMENT, 
            layupNames=('plate_0.partition_square', ),
            layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True, 
            outputAtPlyBottom=False, rebar=EXCLUDE)
        
        self.model.FieldOutputRequest(name='Layup-Output-1', 
            createStepName='Loading', variables=variables,
            frequency=LAST_INCREMENT, 
            layupNames=('plate_1.partition_square', ),
            layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True, 
            outputAtPlyBottom=False, rebar=EXCLUDE)
        
        self.model.FieldOutputRequest(name='Layup-Output-remainder-0', 
            createStepName='Loading', variables=variables,
            frequency=LAST_INCREMENT, 
            layupNames=('plate_0.remainder', ),
            layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True, 
            outputAtPlyBottom=False, rebar=EXCLUDE)
        
        self.model.FieldOutputRequest(name='Layup-Output-remainder-1', 
            createStepName='Loading', variables=variables,
            frequency=LAST_INCREMENT, 
            layupNames=('plate_1.remainder', ),
            layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True, 
            outputAtPlyBottom=False, rebar=EXCLUDE)
        
        self.model.HistoryOutputRequest(name='H-Fastener-0', 
            createStepName='Loading', variables=('CTF1', 'CTF2', 'CTF3', 'CTM1', 
            'CTM2', 'CTM3', 'CU1', 'CU2', 'CU3', 'CUR1', 'CUR2', 'CUR3'), 
            frequency=LAST_INCREMENT, fasteners='Fastener-0', sectionPoints=DEFAULT, 
            rebar=EXCLUDE)
        
        self.model.FieldOutputRequest(name='Fastener-Output-0', 
            createStepName='Loading', variables=('CTF', 'CU'), 
            frequency=LAST_INCREMENT, fasteners='Fastener-0', sectionPoints=DEFAULT, 
            rebar=EXCLUDE)
        
    def setup_jobs(self):
        '''
        Define jobs
        '''
        mdb.Job(name=self.name_job, model=str(self.name_model), description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, 
            memory=self.pRun['memory_max_percentage'], 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, 
            numThreadsPerMpiProcess=self.pRun['numThreadsPerMpiProcess'], 
            multiprocessingMode=DEFAULT, 
            numCpus=self.pRun['numCpus'], numDomains=self.pRun['numCpus'], numGPUs=0)


def extract_field_S4R(name_job, fname_save='specimen-field-S4R.dat'):
    '''
    Extract field output at integration points for shell elements.
    '''
    N_SET = 2
    NAME_INSTANCES = ['PLATE_0', 'PLATE_1']
    NAME_SETS = ['PARTITION_SQUARE', 'PARTITION_SQUARE']
    
    odb = OdbOperation(name_job)

    with open(fname_save, 'w') as f:
        
        f.write('Variables= X Y Z thickness index S11 S22 S12 TSHR13 TSHR23 E11 E22 E12\n')
        
        for i_set in range(N_SET):
        
            name_instance = NAME_INSTANCES[i_set]
        
            element_labels, indices_fieldOutput = odb.get_element_labels_and_indices(
                name_instance=name_instance, name_set=NAME_SETS[i_set])
            coordinates = odb.probe_element_center_coordinate(
                name_instance=name_instance, element_label=element_labels)
        
            f.write('zone T=" %s %s "\n'%(name_instance, NAME_SETS[i_set]))
            
            for i_elem, element_label in enumerate(element_labels):
                
                values_S11 = odb.probe_shell_element_thickness_values(variable='S', component='S11',
                                        name_instance=name_instance, element_label=element_label)
                values_S22 = odb.probe_shell_element_thickness_values(variable='S', component='S22',
                                        name_instance=name_instance, element_label=element_label)
                values_S12 = odb.probe_shell_element_thickness_values(variable='S', component='S12',
                                        name_instance=name_instance, element_label=element_label)
                values_TSHR13 = odb.probe_shell_element_thickness_values(variable='TSHR13', component=None,
                                        name_instance=name_instance, element_label=element_label)
                values_TSHR23 = odb.probe_shell_element_thickness_values(variable='TSHR23', component=None,
                                        name_instance=name_instance, element_label=element_label)
                values_E11 = odb.probe_shell_element_thickness_values(variable='E', component='E11',
                                        name_instance=name_instance, element_label=element_label)
                values_E22 = odb.probe_shell_element_thickness_values(variable='E', component='E22',
                                        name_instance=name_instance, element_label=element_label)
                values_E12 = odb.probe_shell_element_thickness_values(variable='E', component='E12',
                                        name_instance=name_instance, element_label=element_label)

                thickness_distribution = values_S11[:, 0]
                n_thickness = len(thickness_distribution)
            
                for i_thick in range(1, n_thickness-1):
                    for j in range(3):
                        f.write(' %14.6E'%(coordinates[i_elem][j]))
                    f.write(' %14.6E'%(thickness_distribution[i_thick]))
                    f.write(' %d'%(indices_fieldOutput[i_elem]))
                    f.write(' %14.6E'%(values_S11[i_thick,1]))
                    f.write(' %14.6E'%(values_S22[i_thick,1]))
                    f.write(' %14.6E'%(values_S12[i_thick,1]))
                    f.write(' %14.6E'%(values_TSHR13[i_thick-1,1]))
                    f.write(' %14.6E'%(values_TSHR23[i_thick-1,1]))
                    f.write(' %14.6E'%(values_E11[i_thick,1]))
                    f.write(' %14.6E'%(values_E22[i_thick,1]))
                    f.write(' %14.6E'%(values_E12[i_thick,1]))
                    f.write('\n')

def extract_mid_plane_strain(name_job, fname_save='specimen-mid-plane-strain-S4R.dat'):
    
    N_SET = 2
    NAME_INSTANCES = ['PLATE_0', 'PLATE_1']
    NAME_SETS = ['PARTITION_SQUARE', 'PARTITION_SQUARE']
    
    odb = OdbOperation(name_job)
    
    with open(fname_save, 'w') as f:
        
        f.write('Variables= X Y Z index N11 N22 N12 M11 M22 M12')
        f.write(' epsilon11 epsilon22 epsilon12 kappa11 kappa22 kappa12\n')
    
        for i_set in range(N_SET):
            
            name_instance = NAME_INSTANCES[i_set]
            
            element_labels, indices_fieldOutput = odb.get_element_labels_and_indices(
                name_instance=name_instance, name_set=NAME_SETS[i_set])
            coordinates = odb.probe_element_center_coordinate(
                name_instance=name_instance, element_label=element_labels)

            _element_labels, value_SE = odb.probe_element_set_values(
                step='Loading', frame=-1, variable='SE', component=None,
                name_instance=name_instance, name_set=NAME_SETS[i_set])
            
            _, value_SF = odb.probe_element_set_values(
                step='Loading', frame=-1, variable='SF', component=None,
                name_instance=name_instance, name_set=NAME_SETS[i_set])

            _, value_SM = odb.probe_element_set_values(
                step='Loading', frame=-1, variable='SM', component=None,
                name_instance=name_instance, name_set=NAME_SETS[i_set])
            
            _, value_SK = odb.probe_element_set_values(
                step='Loading', frame=-1, variable='SK', component=None,
                name_instance=name_instance, name_set=NAME_SETS[i_set])
            
            n_elements = len(element_labels)
            f.write('zone T=" %s %s " I= %d\n'%(name_instance, NAME_SETS[i_set], n_elements))

            for i_elem in range(n_elements):
                if _element_labels[i_elem] != element_labels[i_elem]:
                    raise ValueError('Element labels do not match: %d != %d'%(
                        _element_labels[i_elem], element_labels[i_elem]))
                for j in range(3):
                    f.write(' %14.6E'%(coordinates[i_elem][j]))
                f.write(' %d'%(indices_fieldOutput[i_elem]))
                f.write(' %14.6E'%(value_SF[i_elem, 0])) # SF1 = N11
                f.write(' %14.6E'%(value_SF[i_elem, 1])) # SF2 = N22
                f.write(' %14.6E'%(value_SF[i_elem, 3])) # SF6 = N12
                f.write(' %14.6E'%(value_SM[i_elem, 1])) # SM1 = M11
                f.write(' %14.6E'%(value_SM[i_elem, 0])) # SM2 = M22
                f.write(' %14.6E'%(value_SM[i_elem, 2])) # SM3 = M12
                f.write(' %14.6E'%(value_SE[i_elem, 0])) # SE1 = epsilon11
                f.write(' %14.6E'%(value_SE[i_elem, 1])) # SE2 = epsilon22
                f.write(' %14.6E'%(value_SE[i_elem, 3])) # SE6 = epsilon12
                f.write(' %14.6E'%(value_SK[i_elem, 1])) # SK1 = kappa11
                f.write(' %14.6E'%(value_SK[i_elem, 0])) # SK2 = kappa22
                f.write(' %14.6E'%(value_SK[i_elem, 2])) # SK3 = kappa12
                f.write('\n')


if __name__ == '__main__':

    if os.path.exists('parameters.json'):
        print('>>> Found [parameters.json], use it for the model setup.')
        with open('parameters.json', 'r') as f:
            parameters = json.load(f)
    else:
        with open('default-parameters.json', 'r') as f:
            parameters = json.load(f)

    pGeo = parameters['pGeo']
    pMesh = parameters['pMesh']
    pRun = parameters['pRun']
    
    index_run = parameters['index_run']
    index_case = parameters['index_case']
    displacement = parameters['displacement']
    
    #* Implicit modelling dict for mesh-independent fasteners
    pMesh['ImplicitModelling'] = {
        "width_partition": pGeo['fasteners'][0]['r_washer']*1.5,
        "vf_hole": 0.4,
        "E11": pMesh['E11'],
        "E22": pMesh['E22'],
        "E33": pMesh['E33'],
        "G12": pMesh['G12'],
        "G13": pMesh['G13'],
        "G23": pMesh['G23'],
        "nu12": pMesh['nu12'],
        "nu13": pMesh['nu13'],
        "nu23": pMesh['nu23']
    }

    #* Build model

    name_job = 'Job_MIF_%d_%d'%(index_run, index_case)

    model = SingleLapBoltedJoint(name_job, pGeo, pMesh, pRun, displacement=displacement)
    model.build()
    model.set_view()
    model.save_cae('MIF_%d_%d.cae'%(index_run, index_case))
    
    if not parameters['not_run_job']:
        
        print('>>> Running job [%s]...'%(name_job))
        model.write_job_inp(model.name_job)
        model.submit_job(model.name_job, only_data_check=False)

        #* Post process
        print('>>> Extracting results from job [%s]...'%(name_job))
        odb = OdbOperation(model.name_job)

        with open(model.name_job+'-RF.dat', 'w') as f:
            
            _, u_RP = odb.probe_node_set_values(step='Loading', frame=-1, variable='U',
                        component=None, name_instance=None, name_set='RF_LOAD')
            u_RP = u_RP[0] # only one node in the set

            _, rf_RP = odb.probe_node_set_values(step='Loading', frame=-1, variable='RF',
                        component=None, name_instance=None, name_set='RF_LOAD')
            rf_RP = rf_RP[0] # only one node in the set

            for i in range(3):
                f.write('%s_RF%d  %20.6E \n'%(model.label_rp, i+1, rf_RP[i]))
            for i in range(3):
                f.write('%s_U%d   %20.6E \n'%(model.label_rp, i+1, u_RP[i]))

            element_labels, ctf_0 = odb.probe_element_set_values(step='Loading', frame=-1, variable='CTF',
                        component=None, name_instance=None, name_set='_FASTENER-0_PF_',
                        position=WHOLE_ELEMENT) # Position is WHOLE_ELEMENT for CTF, CTM, CU
            ctf_0 = ctf_0[0] # only one element in the set
            for i in range(3):
                f.write('FASTENER%d_CTF%d   %20.6E \n'%(0, i+1, ctf_0[i]))

            element_labels, ctm_0 = odb.probe_element_set_values(step='Loading', frame=-1, variable='CTM',
                        component=None, name_instance=None, name_set='_FASTENER-0_PF_',
                        position=WHOLE_ELEMENT) # Position is WHOLE_ELEMENT for CTF, CTM, CU
            ctm_0 = ctm_0[0] # only one element in the set
            for i in range(3):
                f.write('FASTENER%d_CTM%d   %20.6E \n'%(0, i+1, ctm_0[i]))

        extract_field_S4R(name_job=name_job, fname_save=name_job+'-field.dat')
        extract_mid_plane_strain(name_job=name_job, fname_save=name_job+'-mid-plane.dat')
            
