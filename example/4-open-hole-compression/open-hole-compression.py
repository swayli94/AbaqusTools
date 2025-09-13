'''
Open Hole Compression Test of a Composite Plate using Abaqus/CAE Python Scripting

- The plate is compressed in the y-axis direction.
- The thickness direction is the z-axis.
- The plate is clamped at the zx-plane (bottom face in the y-axis direction).
- The displacement is applied at the top surface in the y-axis direction.

'''
import os
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


class Plate(Part):
    '''
    Plate with an open hole.
    '''
    def __init__(self, model, pGeo, pMesh):
        
        super(Plate, self).__init__(model, pGeo, pMesh)
        
        self.name_part = 'plate'

        #* Attributes
        self.len_x = self.pGeo['len_x_plate']
        self.len_y = self.pGeo['len_y_plate']
        self.thk_z = self.pGeo['thk_z_plate']
        
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
        
        if self.pMesh['ply_model'] == 'PlyByPly':
            
            self.create_partition_ply()
            self.loop_over_plies()
        
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
            mySkt.Line(point1=tuple(sketch_points[i-1,:]), point2=tuple(sketch_points[i,:]))

        o0 = np.array([self.xc_hole, self.yc_hole])
        dd = np.array([self.r_hole, self.r_hole]) * np.sqrt(2)*0.5
        mySkt.CircleByCenterPerimeter(center=tuple(o0), point1=tuple(o0+dd))
        
        #* Create the partition circle-square sketch in x-y plane
        mySkt = self.model.ConstrainedSketch(name='partition_top_view', sheetSize=200)
        
        o0 = np.array([self.xc_hole, self.yc_hole])
        dd = np.array([self.r_partition, self.r_partition]) * np.sqrt(2)*0.5
        mySkt.CircleByCenterPerimeter(center=tuple(o0), point1=tuple(o0+dd))

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
        
        self.create_datum_csys_3p(myPrt, 'csys_plate', origin=[0.0, 0.0, 0.0], dx=[1, 0, 0], dy=[0, 1, 0])
        
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
        myPrt.BaseSolidExtrude(sketch=mySkt, depth=self.thk_z)

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

        faces = self.get_faces(myPrt, (0, 0.5*self.len_y, 0.5*self.thk_z))
        myPrt.Surface(side1Faces=faces, name='face_x0')
        
        faces = self.get_faces(myPrt, (self.len_x, 0.5*self.len_y, 0.5*self.thk_z))
        myPrt.Surface(side1Faces=faces, name='face_x1')

        faces = self.get_faces(myPrt, (0.5*self.len_x, 0.0, 0.5*self.thk_z))
        myPrt.Surface(side1Faces=faces, name='face_y0')
        
        faces = self.get_faces(myPrt, (0.5*self.len_x, self.len_y, 0.5*self.thk_z))
        myPrt.Surface(side1Faces=faces, name='face_y1')
        
        faces = self.get_faces(myPrt, (pt_x, pt_y, 0.0))
        myPrt.Surface(side1Faces=faces, name='face_z0')
        
        faces = self.get_faces(myPrt, (pt_x, pt_y, self.thk_z))
        myPrt.Surface(side1Faces=faces, name='face_z1')
        
        faces = self.get_faces(myPrt, (self.xc_hole + self.r_hole, self.yc_hole, 0.5*self.thk_z))
        myPrt.Surface(side1Faces=faces, name='face_hole')

    def create_set(self):
        
        pt_x = 0.5*(self.xc_hole - self.r_hole)
        pt_y = 0.5*(self.yc_hole - self.r_hole)
        
        myPrt = self.model.parts[self.name_part]
        myPrt.Set(cells=myPrt.cells, name='all') 

        faces = self.get_faces(myPrt, (0.5*self.len_x, 0.0, 0.5*self.thk_z))
        myPrt.Set(faces=faces, name='face_y0')
        
        faces = self.get_faces(myPrt, (0.5*self.len_x, self.len_y, 0.5*self.thk_z))
        myPrt.Set(faces=faces, name='face_y1')

        faces = self.get_faces(myPrt, (pt_x, pt_y, 0.0))
        myPrt.Set(faces=faces, name='face_z0')

        faces = self.get_faces(myPrt, (pt_x, pt_y, self.thk_z))
        myPrt.Set(faces=faces, name='face_z1')
        
        faces = self.get_faces(myPrt, (self.xc_hole + self.r_hole, self.yc_hole, 0.5*self.thk_z))
        myPrt.Set(faces=faces, name='face_hole')
        
        edges = self.get_edges(myPrt, (0.0, 0.0, 0.5*self.thk_z))
        myPrt.Set(edges=edges, name='edge_z_x0y0')
        
        edges = self.get_edges(myPrt, (self.xc_hole + self.r_hole, self.yc_hole, 0.0))
        myPrt.Set(edges=edges, name='edge_hole_z0')
        
        edges = self.get_edges(myPrt, (self.xc_hole + self.r_hole, self.yc_hole, self.thk_z))
        myPrt.Set(edges=edges, name='edge_hole_z1')
    
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
        cells = self.get_cells(myPrt, (self.xc_hole + dd, self.yc_hole, 0.5*self.thk_z))
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
        cells = self.get_cells(myPrt, (self.xc_hole + dd, self.yc_hole, 0.5*self.thk_z))
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
    
    #* Partition, surface, set for ply-by-ply modeling
    
    def create_partition_ply(self):
        '''
        Create partition for each ply
        '''
        myPrt = self.model.parts[self.name_part]
        
        num_ply = self.get_num_ply()
        
        z_top = self.thk_z
        z_bottom = 0.0
        
        for i in range(num_ply-1):
            
            r = (i+1.0)/num_ply
            z = (1-r)*z_bottom + r*z_top
            
            myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.cells, 
                point1=tuple([0, 0, z]), point2=tuple([1, 0, z]), point3=tuple([0, 1, z]))
        
    def loop_over_plies(self):
        '''
        Loop over plies: seed edge, and create sets
        '''
        myPrt = self.model.parts[self.name_part]
        
        #* Ply parameters
        num_ply = self.get_num_ply()
        
        #* Stack direction of plate,
        #* the reference face is the top surface, the stacking direction is from bottom to top,
        #* The 1st ply is in the bottom surface (z0)
        z_top = self.thk_z
        z_bottom = 0.0
        
        t0 = time.time()
        for i_ply in range(num_ply):
            
            t1 = time.time()
            r0 = (i_ply*1.0)/num_ply
            r1 = (i_ply+1.0)/num_ply
            z0 = (1-r0)*z_bottom + r0*z_top
            z1 = (1-r1)*z_bottom + r1*z_top

            self._seed_edge_ply(myPrt, z0, z1)
            
            self._create_set_ply(z0, z1, i_ply)
            
            t2 = time.time()
            print('>>> Seeding ply %2d of [%s], t= %.1f s'%(i_ply+1, self.name_part, t2-t1))
            
        print('>>> Seeding [%s], t= %.1f min'%(self.name_part, (t2-t0)/60.0))

    #* Meshing
    
    def set_seeding(self):

        myPrt = self.model.parts[self.name_part]
        myPrt.seedPart(size=self.pMesh['plate_seedPart_size'], 
                        deviationFactor=0.1, minSizeFactor=0.1)
        
        if self.pMesh['ply_model'] == 'ContinuumShell':
        
            self._seed_edge_cs()
        
        elif self.pMesh['ply_model'] == 'PlyByPly':
            
            #* Seeding in the loop_over_plies
            pass
        
        else:
            
            raise ValueError('Unknown ply model: %s'%(self.pMesh['ply_model']))

    def create_mesh(self):
        
        if self.pMesh['ply_model'] == 'ContinuumShell':
        
            self._create_mesh_cs()
        
        elif self.pMesh['ply_model'] == 'PlyByPly':
            
            self._create_mesh_ply()
        
        else:
            
            raise ValueError('Unknown ply model: %s'%(self.pMesh['ply_model']))

    def set_element_type(self):
        
        if self.pMesh['ply_model'] == 'ContinuumShell':
        
            self._set_element_type_cs()
        
        elif self.pMesh['ply_model'] == 'PlyByPly':
            
            self._set_element_type_ply()
        
        else:
            
            raise ValueError('Unknown ply model: %s'%(self.pMesh['ply_model']))
    
    def set_section_assignment(self):
        
        if self.pMesh['ply_model'] == 'ContinuumShell':
        
            self._set_section_assignment_cs()
        
        elif self.pMesh['ply_model'] == 'PlyByPly':
            
            self._set_section_assignment_ply()
        
        else:
            
            raise ValueError('Unknown ply model: %s'%(self.pMesh['ply_model']))
    
    def get_num_ply(self):
        '''
        Get number of plies in the plate
        
        Returns
        ------------
        num_ply: int
            number of plies
        '''
        return len(self.pMesh['plate_CompositePly_orientationValue']) * (1 + self.pMesh['plate_CompositeLayup_symmetric'])
    
    def get_angle_ply(self, i_ply):
        '''
        Get the angle of ply in the plate
        
        Parameters
        ------------
        i_ply: int
            index of the ply, starts from 0
            
        Returns
        ------------
        angle: float
            the composite ply's orientation angle (degree)
        '''
        layup = self.pMesh['plate_CompositePly_orientationValue']
        
        ii = i_ply
        
        if self.pMesh['plate_CompositeLayup_symmetric'] and i_ply>=len(layup):
        
            ii = 2*len(layup) - i_ply - 1

        return layup[ii]
    
    #* Continuum shell modeling
    
    def _seed_edge_cs(self):
        
        myPrt = self.model.parts[self.name_part]
        
        myPrt.seedEdgeByNumber(edges=myPrt.sets['edge_z_x0y0'].edges, 
                                number=self.pMesh['num_element_thickness'], constraint=FIXED)
        
        self._seed_edge_face_hole_radial(myPrt, 0.0, reverse=False)
        self._seed_edge_face_hole_radial(myPrt, self.thk_z, reverse=True)
        
        self._seed_edge_face_circumferential_partition(myPrt, 0.0)
        self._seed_edge_face_circumferential_partition(myPrt, self.thk_z)

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
                    
    def _create_mesh_cs(self):

        #* Stack direction of plate,
        #* the reference face is the top surface, the stacking direction is from bottom to top,
        #* The 1st ply is in the bottom surface (z0)
        #* The top face is the z1 face

        myPrt = self.model.parts[self.name_part]
        # myPrt.setMeshControls(regions=myPrt.cells, technique=SWEEP, algorithm=ADVANCING_FRONT)
        myPrt.setMeshControls(regions=myPrt.cells, elemShape=HEX)
        myPrt.assignStackDirection(referenceRegion=myPrt.surfaces['face_z1'].faces[0], cells=myPrt.cells)
        myPrt.generateMesh()
    
    def _set_element_type_cs(self):

        myPrt = self.model.parts[self.name_part]
        self.set_element_type_of_part(myPrt, kind='continuum shell')
    
    def _set_section_assignment_cs(self):
        
        myPrt = self.model.parts[self.name_part]
        
        self.set_CompositeLayup_of_set(myPrt, 
                name_set=       'all', 
                total_thickness=self.thk_z, 
                ply_angle=      self.pMesh['plate_CompositePly_orientationValue'],
                eNum_thickness= self.pMesh['num_element_thickness'],
                symmetric=      self.pMesh['plate_CompositeLayup_symmetric'],
                numIntPoints=   self.pMesh['plate_CompositePly_numIntPoints'],
                name_csys_datum='csys_plate')
                
    #* Ply-by-ply modeling 

    def _seed_edge_ply(self, myPrt, z0, z1):
        '''
        Seed edges on one face of the ply partition
        
        Parameters
        ------------------
        myPrt: Abaqus part object
            part of the plate

        z0, z1: float
            z coordinates of the ply partition faces
        '''
        #* Thickness direction edge (edge_z_x0y1)
        edges = self.get_edges(myPrt, (0.0, 0.0, 0.5*(z0+z1)))
        myPrt.seedEdgeByNumber(edges=edges, number=self.pMesh['num_element_thickness'], constraint=FIXED)

        #* Face edges
        self._seed_edge_face_hole_radial(myPrt, z0, reverse=False)
        self._seed_edge_face_circumferential_partition(myPrt, z0)
 
        if z1 == self.thk_z:
            
            self._seed_edge_face_hole_radial(myPrt, z1, reverse=True)
            self._seed_edge_face_circumferential_partition(myPrt, z1)

    def _create_set_ply(self, z0, z1, i_ply):
        '''
        Create set for each ply
        
        Parameters
        ------------------
        z0, z1: float
            z coordinates of the ply partition faces

        i_ply: int
            index of the ply
        '''
        z_mid = 0.5*(z0+z1)
        width = 0.5*self.width_partition

        EPSILON = 0.001
        ANGLE_INCREMENT= 0.5*np.pi

        #* Cells around the hole
        dc = 0.5*(self.r_hole + self.r_partition)
        ds = 0.5*(self.r_partition + self.width_partition*0.5)

        points=[]
        for i in range (4):
            
            angle= ANGLE_INCREMENT * i
            
            x = self.xc_hole + dc*np.sin(angle)
            y = self.yc_hole + dc*np.cos(angle)
            points.append((x,y,z_mid))
            
            x = self.xc_hole + ds*np.sin(angle)
            y = self.yc_hole + ds*np.cos(angle)
            points.append((x,y,z_mid))
            
        #* Cells of rectangular blocks
        points += [
            (self.xc_hole - width - EPSILON,    self.yc_hole - width - EPSILON, z_mid),
            (self.xc_hole,                      self.yc_hole - width - EPSILON, z_mid),
            (self.xc_hole + width + EPSILON,    self.yc_hole - width - EPSILON, z_mid),
            (self.xc_hole - width - EPSILON,    self.yc_hole,                   z_mid),
            (self.xc_hole + width + EPSILON,    self.yc_hole,                   z_mid),
            (self.xc_hole - width - EPSILON,    self.yc_hole + width + EPSILON, z_mid),
            (self.xc_hole,                      self.yc_hole + width + EPSILON, z_mid),
            (self.xc_hole + width + EPSILON,    self.yc_hole + width + EPSILON, z_mid),
        ]

        #* Create set
        self.create_geometry_set('ply-%d'%(i_ply+1), self.name_part, points, geometry='cell')

    def _create_mesh_ply(self):

        #* Stack direction of plate,
        #* the reference face is the top surface, the stacking direction is from bottom to top,
        #* The 1st ply is in the bottom surface (z0)
        #* The top face is the z1 face

        t0 = time.time()
        
        myPrt = self.model.parts[self.name_part]
        myPrt.setMeshControls(regions=myPrt.cells, elemShape=HEX)
        myPrt.assignStackDirection(referenceRegion=myPrt.surfaces['face_z1'].faces[0], cells=myPrt.cells)
        myPrt.generateMesh()

        t1 = time.time()
        print('>>> Meshing of [%s], t= %.1f min'%(self.name_part, (t1-t0)/60.0))
        
    def _set_element_type_ply(self):
        
        myPrt = self.model.parts[self.name_part]
        
        self.set_element_type_of_part(myPrt, kind='3D stress', hourglassControl='enhanced')
        
    def _set_section_assignment_ply(self):
            
        #* Stack direction of plate,
        #* the reference face is the top surface, the stacking direction is from bottom to top,
        #* The 1st ply is in the bottom surface (z0)
        #* The top face is the z1 face

        myPrt = self.model.parts[self.name_part]
        
        myPrt.SectionAssignment(region=myPrt.sets['all'], sectionName='IM7/8551-7', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
        
        num_ply = self.get_num_ply()

        for i_ply in range(num_ply):
            
            name_set = 'ply-%d'%(i_ply+1)
            
            angle = self.get_angle_ply(i_ply)

            localCsys = self.get_datum_by_name(myPrt, 'csys_plate')
            
            myPrt.MaterialOrientation(region=myPrt.sets[name_set], 
                orientationType=SYSTEM, 
                axis=AXIS_3,                # Additional Rotation Direction
                localCsys=localCsys,        # Orientation by a datum CSYS
                fieldName='', 
                additionalRotationType=ROTATION_ANGLE, 
                additionalRotationField='', 
                angle=angle,                # Additional Rotation angle (degree)
                stackDirection=STACK_3)     # Stacking Direction (STACK_3: bottom to top)


class TestModel(Model):
    
    def __init__(self, name_job, pGeo, pMesh, pRun, displacement=[0.0, -1.0, 0.0]):
        
        super(TestModel,self).__init__(pGeo=pGeo, pMesh=pMesh, pRun=pRun)

        self.name_job = name_job
        self.displacement = displacement
        
    def initialization(self):
        
        self.model = mdb.models[str(self.name_model)]

        if self.pMesh['ply_model'] == 'ContinuumShell':
        
            self.create_material_IM785517(elastic_type='LAMINA')
        
        elif self.pMesh['ply_model'] == 'PlyByPly':
            
            self.create_material_IM785517(elastic_type='ENGINEERING_CONSTANTS')
        
        else:
            
            raise ValueError('Unknown ply model: %s'%(self.pMesh['ply_model']))
        
        self.create_section_IM785517()

    def setup_parts(self):
        
        self.plate = Plate(self.model, self.pGeo, self.pMesh)
        self.plate.build()
    
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
        a = self.rootAssembly
        
        region = a.instances['plate'].sets['face_y0']
        mdb.models['Model-1'].EncastreBC(name='BC-y0', createStepName='Initial', region=region, localCsys=None)
        
        region = a.instances['plate'].sets['face_y1']
        mdb.models['Model-1'].DisplacementBC(name='BC-y1', createStepName='Loading', region=region, 
            u1=self.displacement[0], u2=self.displacement[1], u3=self.displacement[2], 
            ur1=0.0, ur2=0.0, ur3=0.0, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
            fieldName='', localCsys=None)
        
    def setup_outputs(self):
        '''
        Define outputs
        
        https://docs.software.vt.edu/abaqusv2022/English/?show=SIMACAEKERRefMap/simaker-c-fieldoutputrequestpyc.htm
        '''
        #* Every n increments
        frequency = 1
        if self.pRun is not None:
            if 'field_output_frequency' in self.pRun.keys():
                frequency = int(self.pRun['field_output_frequency'])
                if frequency <= 0:
                    frequency = LAST_INCREMENT
                    print('>>> [setup_outputs]:')
                    print('    Only output the last increment')
        
        #* Evenly spaced time intervals
        numIntervals = 0
        if self.pRun is not None:
            if 'field_output_numIntervals' in self.pRun.keys():
                numIntervals = int(self.pRun['field_output_numIntervals'])
                if numIntervals > 0:
                    print('>>> [setup_outputs]:')
                    print('    Output frequency is changed to [Evenly spaced time intervals]')
        
        #* Output variables
        variables1=('S', 'E', 'U', 'RF')

        if self.pMesh['user_subroutine'] == 'UVARM':
            variables1=('S', 'E', 'U', 'RF', 'UVARM')

        if numIntervals <= 0:
            self.model.fieldOutputRequests['F-Output-1'].setValues(variables=variables1, frequency=frequency)
        else:
            self.model.fieldOutputRequests['F-Output-1'].setValues(variables=variables1, numIntervals=numIntervals)
        
        self.setup_field_outputs_composite_layup(frequency=0)
        
    def setup_field_outputs_composite_layup(self, frequency=1, numIntervals=0):
        '''
        Define outputs for composite layup domains, 
        which outputs at all section points of a set (section).
        '''
        if frequency <= 0:
            frequency = LAST_INCREMENT
        
        name_output = 'FO-layup'
        
        self.model.FieldOutputRequest(name=name_output, 
            createStepName='Loading', variables=('S', 'E', 'U', 'SDV'), frequency=frequency,
            layupNames=('plate.partition_circle', ), 
            layupLocationMethod=ALL_LOCATIONS, rebar=EXCLUDE)
        
        if numIntervals > 0:
            self.model.fieldOutputRequests[name_output].setValues(numIntervals=numIntervals)

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

    pGeo = {
        'len_x_plate': 60,
        'len_y_plate': 100,
        'thk_z_plate': 4,
        'xr_hole_center': 0.5,
        'yr_hole_center': 0.5,
        'r_hole': 10,
    }

    pMesh = {
        'ply_model': 'PlyByPly',        # PlyByPly, ContinuumShell
        'failure_model': None,          # 'LaRC05', None
        'user_subroutine': None,        # 'UVARM', None
        'composite_ply_thickness': 0.1,
        'num_element_thickness': 1,
        
        'radius_ratio_partition_circle': 1.5,
        'hole_circumferential_num_seedEdgeByNumber': 16,
        'circle_radial_bias_seedEdgeByBias': 5.0,
        'circle_radial_num_seedEdgeByBias': 16,
        'square_radial_bias_seedEdgeByBias': 3.0,
        'square_radial_num_seedEdgeByBias': 8,
        
        'plate_seedPart_size': 2.0,
        'plate_CompositeLayup_symmetric': True,
        'plate_CompositePly_numIntPoints': 3,
        'plate_CompositePly_orientationValue':[ 45, -45, 0, 0, 45, 90, -45, 0, 0, 90,
                                                45, -45, 0, 0, 45, 90, -45, 0, 0, 90],
    }

    pRun = {
        
        # Static Analysis Step
        
        'timePeriod': 1.0,
        'maxNumInc': 10000000,
        'initialInc': 0.01,
        'minInc': 1e-15,
        'maxInc': 0.1,
        'nlgeom': False,
        
        # Job
        
        'memory_max_percentage': 90,
        'numCpus': 4,
        'numThreadsPerMpiProcess': 1,
        
        # Expert
        
        'contact_friction_coef': 0.5,
        'field_output_frequency': 0,
        'field_output_numIntervals': 0
        
    }

    strain = - 1E-3
    
    #* Load user-input parameters
    if True:
        
        if os.path.exists('input.txt'):
            
            with open('input.txt', 'r') as f:
                
                lines = f.readlines()
                
                pGeo['len_x_plate'] = float(lines[0].split()[1])
                pGeo['len_y_plate'] = float(lines[1].split()[1])
                pGeo['thk_z_plate'] = float(lines[2].split()[1])
                pGeo['xr_hole_center'] = float(lines[3].split()[1])
                pGeo['yr_hole_center'] = float(lines[4].split()[1])
                pGeo['r_hole'] = float(lines[5].split()[1])
                
                strain = float(lines[6].split()[1])
                
                index_layup = int(float(lines[7].split()[1]))
                
            #* Update parameters
            
            pGeo['thk_z_plate'], n_ply = LayupParameters.rounding_thickness(
                                    pGeo['thk_z_plate'], pMesh['composite_ply_thickness'])
            
            pMesh['plate_CompositePly_orientationValue'] = LayupParameters.candidate_composite_layup(
                                    n_ply, index_layup)
       
    #* Build model
    model = TestModel('Job_OHT', pGeo, pMesh, pRun, displacement=[0.0, strain*pGeo['len_y_plate'], 0.0])
    
    model.build()
    model.save_cae('OHT.cae')
    model.write_job_inp()
    
    if pMesh['user_subroutine'] == 'UVARM':
        model.write_IM785517_property_table_inp(method='UVARM', fname_input=model.name_job+'.inp')
