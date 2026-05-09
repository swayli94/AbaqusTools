'''
A laminate with 3D stress elements (S4R).
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


class OpenHolePlate(Part):
    '''
    Plate with an open hole.
    '''
    def __init__(self, model, pGeo, pMesh):
        
        super(OpenHolePlate, self).__init__(model, pGeo, pMesh)
        
        self.name_part = 'plate'

        #* Attributes
        self.length_x = self.pGeo['len_x_plate']
        self.length_y = self.pGeo['len_y_plate']
        self.length_z = self.pGeo['len_z_plate']
        
        self.xc_hole = self.pGeo['xr_hole_center'] * self.length_x
        self.yc_hole = self.pGeo['yr_hole_center'] * self.length_y
        self.r_hole = self.pGeo['r_hole']
        
        self.r_partition, self.width_partition = self._cal_partition_dimensions(
            len_x=self.length_x, len_y=self.length_y, r_hole=self.r_hole,
            xc_hole=self.xc_hole, yc_hole=self.yc_hole,
            radius_ratio_partition_circle=self.pMesh['radius_ratio_partition_circle'])
        
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

    @staticmethod
    def _cal_partition_dimensions(len_x, len_y, r_hole,
                            xc_hole, yc_hole,
                            radius_ratio_partition_circle,
                            ratio_circle=0.3, ratio_square=0.6):
        '''
        Calculate the radius of partition circle around the open hole,
        and the length of partition square around the partition circle.
        '''
        dx0 = xc_hole - r_hole
        dx1 = len_x - xc_hole - r_hole
        dy0 = yc_hole - r_hole
        dy1 = len_y - yc_hole - r_hole
        
        r_partition = min(
            radius_ratio_partition_circle*r_hole, 
            r_hole + ratio_circle*dx0, r_hole + ratio_circle*dx1, 
            r_hole + ratio_circle*dy0, r_hole + ratio_circle*dy1
            )
        
        width_partition = 2* min(
            r_partition*2 - r_hole, 
            r_hole + ratio_square*dx0, r_hole + ratio_square*dx1,
            r_hole + ratio_square*dy0, r_hole + ratio_square*dy1
            )
        
        return r_partition, width_partition

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
            
        o0 = np.array([self.xc_hole, self.yc_hole])
        dd = np.array([self.r_hole, self.r_hole]) * np.sqrt(2)*0.5
        mySkt.CircleByCenterPerimeter(center=tuple(o0),
                                      point1=tuple(o0+dd))
        
        #* Create the partition circle-square sketch in x-y plane
        mySkt = self.model.ConstrainedSketch(name='partition_top_view', sheetSize=200)
        
        o0 = np.array([self.xc_hole, self.yc_hole])
        dd = np.array([self.r_partition, self.r_partition]) * np.sqrt(2)*0.5
        mySkt.CircleByCenterPerimeter(center=tuple(o0), 
                                      point1=tuple(o0+dd))

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
        
        xc_hole = self.xc_hole
        yc_hole = self.yc_hole
        r_hole = self.r_hole
        lx = self.length_x
        ly = self.length_y
        lz = self.length_z
        pt_x = 0.5*(xc_hole - r_hole)
        pt_y = 0.5*(yc_hole - r_hole)
        
        myPrt = self.model.parts[self.name_part]
        myPrt.Set(faces=myPrt.faces, name='all')
        
        self.create_geometry_set('edge_x0', (0.0,    0.5*ly, 0.0), geometry='edge')
        self.create_geometry_set('edge_x1', (lx,     0.5*ly, 0.0), geometry='edge')
        self.create_geometry_set('edge_y0', (0.5*lx, 0.0,    0.0), geometry='edge')
        self.create_geometry_set('edge_y1', (0.5*lx, ly,     0.0), geometry='edge')
        self.create_geometry_set('edge_hole', (xc_hole + r_hole, yc_hole, 0.0), geometry='edge')

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
        
        #* Partition face by sketch
        transform = myPrt.MakeSketchTransform(
            sketchPlane=self.get_datum_by_name(myPrt, 'XYPLANE'),
            sketchUpEdge=self.get_datum_by_name(myPrt, 'XAXIS'), 
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 0.0))
    
        mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=200, transform=transform)
        mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)

        myPrt.projectReferencesOntoSketch(sketch=mySkt, filter=COPLANAR_EDGES)
        
        mySkt.retrieveSketch(sketch=self.model.sketches['partition_top_view'])

        myPrt.PartitionFaceBySketch(
                sketchUpEdge=self.get_datum_by_name(myPrt, 'XAXIS'), 
                faces=myPrt.surfaces['face'].faces,
                sketchOrientation=BOTTOM, sketch=mySkt)

        del self.model.sketches['__profile__']

        dd = 0.5*(self.r_partition+self.r_hole)
        self.create_geometry_set('partition_circle',
                    (self.xc_hole+dd, self.yc_hole, 0.0), geometry='face')

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

        dd = 0.5*(self.width_partition*0.5+self.r_partition)
        self.create_geometry_set('partition_square',
                    (self.xc_hole+dd, self.yc_hole, 0.0), geometry='face')
        
        #* Partition cell by diagonal planes
        myPrt.PartitionFaceByShortestPath(myPrt.sets['partition_circle'].faces, 
            point1=(x0, y0, 0.0), point2=(x1, y1, 0.0))
        myPrt.PartitionFaceByShortestPath(myPrt.sets['partition_square'].faces, 
            point1=(x0, y0, 0.0), point2=(x1, y1, 0.0))
        
        myPrt.PartitionFaceByShortestPath(myPrt.sets['partition_circle'].faces, 
            point1=(x0, y1, 0.0), point2=(x1, y0, 0.0))
        myPrt.PartitionFaceByShortestPath(myPrt.sets['partition_square'].faces, 
            point1=(x0, y1, 0.0), point2=(x1, y0, 0.0))
    
    #* Meshing
    
    def set_seeding(self):

        myPrt = self.model.parts[self.name_part]
        myPrt.seedPart(size=self.pMesh['plate_seedPart_size'], 
                        deviationFactor=0.1, minSizeFactor=0.1)
        
        self._seed_edge_face_hole_radial(myPrt, 0.0, reverse=True)
        
        self._seed_edge_face_circumferential_partition(myPrt, 0.0)
        
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
                name_set=       'all', 
                total_thickness=self.length_z, 
                ply_angle=      self.pMesh['plate_CompositePly_orientationValue'],
                eNum_thickness= self.pMesh['num_element_thickness'],
                symmetric=      self.pMesh['plate_CompositeLayup_symmetric'],
                numIntPoints=   self.pMesh['plate_CompositePly_numIntPoints'],
                name_csys_datum='csys_plate',
                material='orthotropic',
                elementType='shell')
        
    def _seed_edge_face_hole_radial(self, myPrt, z, reverse=False):
        '''
        Seed the edges around the hole in radial direction in one face.
        '''
        ratio_circle  = self.pMesh['circle_radial_bias_seedEdgeByBias']
        number_circle = self.pMesh['circle_radial_num_seedEdgeByBias']
        ratio_square  = self.pMesh['square_radial_bias_seedEdgeByBias']
        number_square = self.pMesh['square_radial_num_seedEdgeByBias']

        ANGLE_INCREMENT= 0.25*np.pi

        dc = 0.5*(self.r_hole + self.r_partition)
        ds = 0.5*(self.r_partition + self.width_partition*0.5)

        for i in range (4):
            
            angle= ANGLE_INCREMENT * (2*i+1)
            
            x_c = self.xc_hole + dc*np.sin(angle)
            y_c = self.yc_hole + dc*np.cos(angle)

            x_s = self.xc_hole + ds*np.sin(angle)
            y_s = self.yc_hole + ds*np.cos(angle)

            edges_c = self.get_edges(myPrt, (x_c,y_c,z), getClosest=False)
            edges_s = self.get_edges(myPrt, (x_s,y_s,z), getClosest=False)
            
            if not reverse:
                
                if i in [0, 1]:
                    myPrt.seedEdgeByBias(biasMethod=SINGLE, end1Edges=edges_c, ratio=ratio_circle, number=number_circle, constraint=FIXED)
                    myPrt.seedEdgeByBias(biasMethod=SINGLE, end1Edges=edges_s, ratio=ratio_square, number=number_square, constraint=FIXED)
                else:
                    myPrt.seedEdgeByBias(biasMethod=SINGLE, end2Edges=edges_c, ratio=ratio_circle, number=number_circle, constraint=FIXED)
                    myPrt.seedEdgeByBias(biasMethod=SINGLE, end2Edges=edges_s, ratio=ratio_square, number=number_square, constraint=FIXED)

            else:

                if i in [0, 1]:
                    myPrt.seedEdgeByBias(biasMethod=SINGLE, end2Edges=edges_c, ratio=ratio_circle, number=number_circle, constraint=FIXED)
                    myPrt.seedEdgeByBias(biasMethod=SINGLE, end2Edges=edges_s, ratio=ratio_square, number=number_square, constraint=FIXED)
                else:
                    myPrt.seedEdgeByBias(biasMethod=SINGLE, end1Edges=edges_c, ratio=ratio_circle, number=number_circle, constraint=FIXED)
                    myPrt.seedEdgeByBias(biasMethod=SINGLE, end1Edges=edges_s, ratio=ratio_square, number=number_square, constraint=FIXED)

    def _seed_edge_face_circumferential_partition(self, myPrt, z):
        '''
        Seed the circumferential edges around the hole in one face.
        '''
        num_circum = self.pMesh['hole_circumferential_num_seedEdgeByNumber']
        
        #* Hole edges
        edges = self.get_edges(myPrt, (self.xc_hole - self.r_hole, self.yc_hole, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        edges = self.get_edges(myPrt, (self.xc_hole + self.r_hole, self.yc_hole, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        edges = self.get_edges(myPrt, (self.xc_hole, self.yc_hole - self.r_hole, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        edges = self.get_edges(myPrt, (self.xc_hole, self.yc_hole + self.r_hole, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        #* Partition circle edges
        edges = self.get_edges(myPrt, (self.xc_hole - self.r_partition, self.yc_hole, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        edges = self.get_edges(myPrt, (self.xc_hole + self.r_partition, self.yc_hole, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        edges = self.get_edges(myPrt, (self.xc_hole, self.yc_hole - self.r_partition, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        edges = self.get_edges(myPrt, (self.xc_hole, self.yc_hole + self.r_partition, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        #* Partition square edges
        edges = self.get_edges(myPrt, (self.xc_hole - 0.5*self.width_partition, self.yc_hole, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 

        edges = self.get_edges(myPrt, (self.xc_hole + 0.5*self.width_partition, self.yc_hole, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        edges = self.get_edges(myPrt, (self.xc_hole, self.yc_hole - 0.5*self.width_partition, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
        
        edges = self.get_edges(myPrt, (self.xc_hole, self.yc_hole + 0.5*self.width_partition, z))
        myPrt.seedEdgeByNumber(edges=edges, number=num_circum, constraint=FIXED) 
       

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

    def setup_parts(self):
        
        self.plate = OpenHolePlate(self.model, self.pGeo, self.pMesh)
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
        
        variables = ('S', 'E')
        if 'failure_model' in self.pMesh:
            if self.pMesh['failure_model'] == 'Hashin':
                variables_hashin = ('DMICRT', 'HSNFTCRT', 'HSNFCCRT', 'HSNMTCRT', 'HSNMCCRT')
                variables += variables_hashin

        self.model.FieldOutputRequest(name='Layup-Output', 
            createStepName='Loading', variables=variables,
            frequency=LAST_INCREMENT, 
            layupNames=('plate.all', ), 
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


if __name__ == '__main__':

    with open('default-parameters.json', 'r') as f:
        default_parameters = json.load(f)

    pGeo = default_parameters['pGeo']
    pMesh = default_parameters['pMesh']
    pRun = default_parameters['pRun']
    
    name_job = 'Demo_OHP_laminate_UD_C3D8R'

    model = LaminateModel(name_job, pGeo, pMesh, pRun)
    model.build()
    model.set_view()
