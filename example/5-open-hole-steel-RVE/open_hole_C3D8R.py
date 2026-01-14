'''
Steel plate (in x-y plane) with an open hole in the center.
'''
import json
import time
import numpy as np
from AbaqusTools import Part, IS_ABAQUS

if IS_ABAQUS:
    from abaqus import *
    from abaqusConstants import *
    from sketch import *
    import mesh
    
from AbaqusTools import Part, Model
from AbaqusTools.functions import LayupParameters


class OpenHolePlate(Part):
    '''
    Plate (in x-y plane) with an open hole in the center.
    '''
    def __init__(self, model, pGeo, pMesh):
        
        super(OpenHolePlate, self).__init__(model, pGeo, pMesh)
        
        self.name_part = 'plate'

        #* Attributes
        self.len_x = self.pGeo['len_x_plate']
        self.len_y = self.pGeo['len_y_plate']
        self.len_z = self.pGeo['len_z_plate']
        
        self.xc_hole = self.pGeo['xr_hole_center'] * self.len_x
        self.yc_hole = self.pGeo['yr_hole_center'] * self.len_y
        self.r_hole = self.pGeo['r_hole']
        
        self._cal_partition_dimensions()

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
        
        #* Abaqus Module: Mesh
        self.set_seeding()
        self.create_mesh()
        self.set_element_type()
        
        #* Abaqus Module: Property
        self.set_section_assignment()
        
        t1 = time.time()

        print('>>> --------------------')
        print('    [Part: %s] build time = %.1f (min)'%(self.name_part, (t1-t0)/60.0))
        print('>>>')

    def _cal_partition_dimensions(self, ratio_circle=0.3, ratio_square=0.6):
        '''
        Calculate the radius of partition circle around the open hole,
        and the length of partition square around the partition circle.
        '''
        dx0 = self.xc_hole - self.r_hole
        dx1 = self.len_x - self.xc_hole - self.r_hole
        dy0 = self.yc_hole - self.r_hole
        dy1 = self.len_y - self.yc_hole - self.r_hole
        
        self.r_partition = min(
            self.pMesh['radius_ratio_partition_circle']*self.r_hole, 
            self.r_hole + ratio_circle*dx0, self.r_hole + ratio_circle*dx1, 
            self.r_hole + ratio_circle*dy0, self.r_hole + ratio_circle*dy1
            )
        
        self.width_partition = 2* min(
            self.r_partition*2 - self.r_hole, 
            self.r_hole + ratio_square*dx0, self.r_hole + ratio_square*dx1,
            self.r_hole + ratio_square*dy0, self.r_hole + ratio_square*dy1
            )

    def create_sketch(self):
        '''
        Create the sketch in X-Y plane.
        '''
        #* Sketch points in X-Y plane (a rectangle)
        sketch_points = np.zeros((4,2))
        sketch_points[0,:] = [0.0,        0.0]
        sketch_points[1,:] = [self.len_x, 0.0]
        sketch_points[2,:] = [self.len_x, self.len_y]
        sketch_points[3,:] = [0.0,        self.len_y]

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

        #* Part by extrusion
        myPrt.BaseSolidExtrude(sketch=mySkt, depth=self.len_z)

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

        faces = self.get_faces(myPrt, (0, 0.5*self.len_y, 0.5*self.len_z))
        myPrt.Surface(side1Faces=faces, name='face_x0')
        
        faces = self.get_faces(myPrt, (self.len_x, 0.5*self.len_y, 0.5*self.len_z))
        myPrt.Surface(side1Faces=faces, name='face_x1')

        faces = self.get_faces(myPrt, (0.5*self.len_x, 0.0, 0.5*self.len_z))
        myPrt.Surface(side1Faces=faces, name='face_y0')
        
        faces = self.get_faces(myPrt, (0.5*self.len_x, self.len_y, 0.5*self.len_z))
        myPrt.Surface(side1Faces=faces, name='face_y1')
        
        faces = self.get_faces(myPrt, (pt_x, pt_y, 0.0))
        myPrt.Surface(side1Faces=faces, name='face_z0')
        
        faces = self.get_faces(myPrt, (pt_x, pt_y, self.len_z))
        myPrt.Surface(side1Faces=faces, name='face_z1')
        
        faces = self.get_faces(myPrt, (self.xc_hole + self.r_hole, self.yc_hole, 0.5*self.len_z))
        myPrt.Surface(side1Faces=faces, name='face_hole')

    def create_set(self):
        
        xc_hole = self.xc_hole
        yc_hole = self.yc_hole
        r_hole = self.r_hole
        lx = self.len_x
        ly = self.len_y
        lz = self.len_z
        pt_x = 0.5*(xc_hole - r_hole)
        pt_y = 0.5*(yc_hole - r_hole)

        myPrt = self.model.parts[self.name_part]
        myPrt.Set(cells=myPrt.cells, name='all') 

        self.create_geometry_set('face_x0', (0.0,    0.5*ly, 0.5*lz), geometry='face')
        self.create_geometry_set('face_x1', (lx,     0.5*ly, 0.5*lz), geometry='face')
        self.create_geometry_set('face_y0', (0.5*lx, 0.0,    0.5*lz), geometry='face')
        self.create_geometry_set('face_y1', (0.5*lx, ly,     0.5*lz), geometry='face')
        self.create_geometry_set('face_z0', (pt_x,   pt_y,   0.0   ), geometry='face')
        self.create_geometry_set('face_z1', (pt_x,   pt_y,   lz    ), geometry='face')
        self.create_geometry_set('face_hole', (xc_hole + r_hole, yc_hole, 0.5*lz), geometry='face')

        self.create_geometry_set('edge_x_y0z0', (0.5*lx, 0.0, 0.0), geometry='edge')
        self.create_geometry_set('edge_x_y1z0', (0.5*lx, ly,  0.0), geometry='edge')
        self.create_geometry_set('edge_x_y0z1', (0.5*lx, 0.0, lz ), geometry='edge')
        self.create_geometry_set('edge_x_y1z1', (0.5*lx, ly,  lz ), geometry='edge')

        self.create_geometry_set('edge_y_z0x0', (0.0, 0.5*ly, 0.0), geometry='edge')
        self.create_geometry_set('edge_y_z1x0', (0.0, 0.5*ly, lz ), geometry='edge')
        self.create_geometry_set('edge_y_z0x1', (lx,  0.5*ly, 0.0), geometry='edge')
        self.create_geometry_set('edge_y_z1x1', (lx,  0.5*ly, lz ), geometry='edge')

        self.create_geometry_set('edge_z_x0y0', (0.0, 0.0, 0.5*lz), geometry='edge')
        self.create_geometry_set('edge_z_x1y0', (lx,  0.0, 0.5*lz), geometry='edge')
        self.create_geometry_set('edge_z_x0y1', (0.0, ly,  0.5*lz), geometry='edge')
        self.create_geometry_set('edge_z_x1y1', (lx,  ly,  0.5*lz), geometry='edge')
        
        self.create_geometry_set('edge_hole_z0', (xc_hole + r_hole, yc_hole, 0.0), geometry='edge')
        self.create_geometry_set('edge_hole_z1', (xc_hole + r_hole, yc_hole, lz ), geometry='edge')

        self.create_geometry_set('vertex_000', (0.0, 0.0, 0.0), geometry='vertex')
        self.create_geometry_set('vertex_100', (lx,  0.0, 0.0), geometry='vertex')
        self.create_geometry_set('vertex_010', (0.0, ly,  0.0), geometry='vertex')
        self.create_geometry_set('vertex_110', (lx,  ly,  0.0), geometry='vertex')
        self.create_geometry_set('vertex_001', (0.0, 0.0, lz ), geometry='vertex')
        self.create_geometry_set('vertex_101', (lx,  0.0, lz ), geometry='vertex')
        self.create_geometry_set('vertex_011', (0.0, ly,  lz ), geometry='vertex')
        self.create_geometry_set('vertex_111', (lx,  ly,  lz ), geometry='vertex')
    
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
                faces=myPrt.surfaces['face_z0'].faces,
                sketchOrientation=BOTTOM, sketch=mySkt)

        del self.model.sketches['__profile__']

        #* Partition cell by the partition circle
        edges = self.get_edges(myPrt, (self.xc_hole - self.r_partition, self.yc_hole, 0.0))
        myPrt.Set(edges=edges, name='edge_partition_circle')
        
        myPrt.PartitionCellByExtrudeEdge(
            line=self.get_datum_by_name(myPrt, 'ZAXIS'), 
            cells=myPrt.cells, edges=edges, sense=FORWARD)
        
        dd = 0.5*(self.r_partition+self.r_hole)
        cells = self.get_cells(myPrt, (self.xc_hole + dd, self.yc_hole, 0.5*self.len_z))
        myPrt.Set(cells=cells, name='partition_circle') 
        
        #* Partition cell to squares by 4 planes
        x0 = self.xc_hole - 0.5*self.width_partition
        x1 = self.xc_hole + 0.5*self.width_partition
        y0 = self.yc_hole - 0.5*self.width_partition
        y1 = self.yc_hole + 0.5*self.width_partition
        
        myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.cells, 
            point1=(x0, y0, 0.0), point2=(x1, y0, 0.0), point3=(x0, y0, 1.0))
        myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.cells, 
            point1=(x0, y1, 0.0), point2=(x1, y1, 0.0), point3=(x0, y1, 1.0))
        myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.cells, 
            point1=(x0, y0, 0.0), point2=(x0, y1, 0.0), point3=(x0, y0, 1.0))
        myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.cells, 
            point1=(x1, y0, 0.0), point2=(x1, y1, 0.0), point3=(x1, y0, 1.0))
        
        dd = 0.5*(self.width_partition*0.5+self.r_partition)
        cells = self.get_cells(myPrt, (self.xc_hole + dd, self.yc_hole, 0.5*self.len_z))
        myPrt.Set(cells=cells, name='partition_square') 
        
        #* Partition cell by diagonal planes
        myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.sets['partition_circle'].cells, 
            point1=(x0, y0, 0.0), point2=(x1, y1, 0.0), point3=(x0, y0, 1.0))
        myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.sets['partition_circle'].cells, 
            point1=(x0, y1, 0.0), point2=(x1, y0, 0.0), point3=(x0, y1, 1.0))
        
        myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.sets['partition_square'].cells, 
            point1=(x0, y0, 0.0), point2=(x1, y1, 0.0), point3=(x0, y0, 1.0))
        myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.sets['partition_square'].cells, 
            point1=(x0, y1, 0.0), point2=(x1, y0, 0.0), point3=(x0, y1, 1.0))
    
    #* Meshing
    
    def set_seeding(self):

        myPrt = self.model.parts[self.name_part]
        myPrt.seedPart(size=self.pMesh['plate_seedPart_size'], 
                        deviationFactor=0.1, minSizeFactor=0.1)
        
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_z_x0y0'].edges, 
                                number=self.pMesh['num_element_thickness'], constraint=FIXED)
        
        self._seed_edge_face_hole_radial(myPrt, 0.0, reverse=False)
        self._seed_edge_face_hole_radial(myPrt, self.len_z, reverse=True)
        
        self._seed_edge_face_circumferential_partition(myPrt, 0.0)
        self._seed_edge_face_circumferential_partition(myPrt, self.len_z)

    def create_mesh(self):
        
        #* Stack direction of plate,
        #* the reference face is the top surface, the stacking direction is from bottom to top,

        myPrt = self.model.parts[self.name_part]
        myPrt.setMeshControls(regions=myPrt.cells, elemShape=HEX)
        myPrt.assignStackDirection(referenceRegion=myPrt.surfaces['face_z1'].faces[0], cells=myPrt.cells)
        myPrt.generateMesh()

    def set_element_type(self):
        
        myPrt = self.model.parts[self.name_part]
        self.set_element_type_of_part(myPrt, kind='3D stress')
    
    def set_section_assignment(self):
        
        myPrt = self.model.parts[self.name_part]
        
        myPrt.SectionAssignment(region=myPrt.sets['all'], sectionName='Steel', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    #* Seeding edges on outer surfaces

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
    

class OpenHolePlateModel(Model):
    
    def __init__(self, name_job, pGeo, pMesh, pRun):
        
        super(OpenHolePlateModel,self).__init__(pGeo=pGeo, pMesh=pMesh, pRun=pRun)

        self.name_job = name_job
        
    def initialization(self):
        
        self.model = mdb.models[str(self.name_model)]

        self.create_material_steel()

        self.create_section_steel()

    def setup_parts(self):
        
        self.plate = OpenHolePlate(self.model, self.pGeo, self.pMesh)
        self.plate.build()
        
        self.length_x = self.plate.len_x
        self.length_y = self.plate.len_y
        self.length_z = self.plate.len_z
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
        '''
        Define loads and boundary conditions
        '''
        raise NotImplementedError('Not implemented yet.')
        
    def setup_outputs(self):

        self.model.fieldOutputRequests['F-Output-1'].setValues(
            variables=('S', 'E', 'U', 'RF'), 
            frequency=1)

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

    with open('parameters.json', 'r') as f:
        parameters = json.load(f)

    pGeo = parameters['pGeo']
    pMesh = parameters['pMesh']
    pRun = parameters['pRun']

    name_job = 'Demo_OHP_steel_C3D8R'
    
    model = OpenHolePlateModel(name_job, pGeo, pMesh, pRun)
    model.build()
    model.set_view()
