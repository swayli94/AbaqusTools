'''
A laminate with 3D stress elements (S4R).
'''
import time
import numpy as np
from AbaqusTools import Part, IS_ABAQUS

if IS_ABAQUS:
    from abaqus import *
    from abaqusConstants import *
    from sketch import *
    
from AbaqusTools import Part, Model


class ImplicitModellingPlate(Part):
    '''
    Plate with an open hole.
    '''
    def __init__(self, model, pGeo, pMesh):
        
        super(ImplicitModellingPlate, self).__init__(model, pGeo, pMesh)
        
        self.name_part = 'plate'

        #* Attributes
        self.length_x = self.pGeo['len_x_plate']
        self.length_y = self.pGeo['len_y_plate']
        self.length_z = self.pGeo['len_z_plate']
        
        self.xc_hole = self.pGeo['xr_hole_center'] * self.length_x
        self.yc_hole = self.pGeo['yr_hole_center'] * self.length_y
        self.r_hole = self.pGeo['r_hole']
        
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
        
        pt_x = 0.5*(self.xc_hole - self.r_hole)
        pt_y = 0.5*(self.yc_hole - self.r_hole)
        
        myPrt = self.model.parts[self.name_part]

        myPrt = self.model.parts[self.name_part]
        faces = self.get_faces(myPrt, (pt_x, pt_y, 0.0))
        myPrt.Surface(side1Faces=faces, name='face')

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
        self.create_geometry_set('partition_square',
                    (self.xc_hole+dd, self.yc_hole, 0.0), geometry='face')
        
        self.create_geometry_set('edge_partition_x0', (x0, 0.5*self.length_y, 0.0), geometry='edge')
        self.create_geometry_set('edge_partition_x1', (x1, 0.5*self.length_y, 0.0), geometry='edge')
        self.create_geometry_set('edge_partition_y0', (0.5*self.length_x, y0, 0.0), geometry='edge')
        self.create_geometry_set('edge_partition_y1', (0.5*self.length_x, y1, 0.0), geometry='edge')
        
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
        
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_partition_x0'].edges, number=1, constraint=FIXED)
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_partition_x1'].edges, number=1, constraint=FIXED)
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_partition_y0'].edges, number=1, constraint=FIXED)
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_partition_y1'].edges, number=1, constraint=FIXED)

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


class LaminateModel(Model):
    
    def __init__(self, name_job, pGeo, pMesh, pRun):
        
        super(LaminateModel,self).__init__(pGeo=pGeo, pMesh=pMesh, pRun=pRun, name_job=name_job)

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
        
        self.plate = ImplicitModellingPlate(self.model, self.pGeo, self.pMesh)
        self.plate.build()
        
        self.length_x = self.plate.length_x
        self.length_y = self.plate.length_y
        self.length_z = self.plate.length_z
        self.volume_box = self.length_x*self.length_y*self.length_z
        self.volume_hole = np.pi*self.plate.r_hole**2*self.length_z
        self.volume = self.volume_box - self.volume_hole
        
    def setup_assembly(self):
        
        a = self.rootAssembly
        p = self.model.parts['plate']
        a.Instance(name='plate', part=p, dependent=ON)
    
    def setup_steps(self):
        '''
        Define a static analysis step
        '''
        self.create_static_step(
            timePeriod= self.pRun['timePeriod'],
            maxNumInc=  self.pRun['maxNumInc'],
            initialInc= self.pRun['initialInc'],
            minInc=     self.pRun['minInc'],
            maxInc=     self.pRun['maxInc'],
            nlgeom=     self.pRun['nlgeom'],
        )
        
    def setup_loads(self):

        raise NotImplementedError('Not implemented yet.')
        
    def setup_outputs(self):

        self.model.FieldOutputRequest(name='F-Output-1', 
            createStepName='Loading', variables=('S', 'E', 'U', 'RF'),
            frequency=LAST_INCREMENT)
        
        variables = ('S', 'TSHR', 'E')
        if 'failure_model' in self.pMesh:
            if self.pMesh['failure_model'] == 'Hashin':
                variables_hashin = ('DMICRT', 'HSNFTCRT', 'HSNFCCRT', 'HSNMTCRT', 'HSNMCCRT')
                variables += variables_hashin

        self.model.FieldOutputRequest(name='Layup-Output-1', 
            createStepName='Loading', variables=variables,
            frequency=LAST_INCREMENT, 
            layupNames=('plate.partition_square', ),
            layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True, 
            outputAtPlyBottom=False, rebar=EXCLUDE)
        
        self.model.FieldOutputRequest(name='Layup-Output-remainder', 
            createStepName='Loading', variables=variables,
            frequency=LAST_INCREMENT, 
            layupNames=('plate.remainder', ),
            layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True, 
            outputAtPlyBottom=False, rebar=EXCLUDE)
        
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

