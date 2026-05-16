'''
Single lap specimen for bolted joints (C3D8R).
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


class OpenHolePlate(Part):
    '''
    Plate with an open hole.
    '''
    def __init__(self, model, name_part, pMesh,
            length_x, length_y, length_z, xc_hole, yc_hole, r_hole):
        
        super(OpenHolePlate, self).__init__(model, pGeo=None, pMesh=pMesh)
        
        self.name_part = name_part

        #* Attributes
        self.length_x = length_x
        self.length_y = length_y
        self.length_z = length_z

        self.xc_hole = xc_hole
        self.yc_hole = yc_hole
        self.r_hole = r_hole
                
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

        #* Part by extrusion
        myPrt.BaseSolidExtrude(sketch=mySkt, depth=self.length_z)

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

        faces = self.get_faces(myPrt, (0, 0.5*self.length_y, 0.5*self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_x0')
        
        faces = self.get_faces(myPrt, (self.length_x, 0.5*self.length_y, 0.5*self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_x1')

        faces = self.get_faces(myPrt, (0.5*self.length_x, 0.0, 0.5*self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_y0')
        
        faces = self.get_faces(myPrt, (0.5*self.length_x, self.length_y, 0.5*self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_y1')
        
        faces = self.get_faces(myPrt, (pt_x, pt_y, 0.0))
        myPrt.Surface(side1Faces=faces, name='face_z0')
        
        faces = self.get_faces(myPrt, (pt_x, pt_y, self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_z1')

        faces = self.get_faces(myPrt, (self.xc_hole + self.r_hole, self.yc_hole, 0.5*self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_hole')

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
        cells = self.get_cells(myPrt, (self.xc_hole + dd, self.yc_hole, 0.5*self.length_z))
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
        cells = self.get_cells(myPrt, (self.xc_hole + dd, self.yc_hole, 0.5*self.length_z))
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
        
        z_top = self.length_z
        z_bottom = 0.0
        
        for i in range(num_ply-1):
            
            r = (i+1.0)/num_ply
            z = (1-r)*z_bottom + r*z_top
            
            myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.cells, 
                point1=tuple([0, 0, z]),
                point2=tuple([1, 0, z]),
                point3=tuple([0, 1, z]))
        
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
        z_top = self.length_z
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
        
    def create_mesh(self):
        
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

    def set_element_type(self):
        
        myPrt = self.model.parts[self.name_part]
        
        self.set_element_type_of_part(myPrt, kind='3D stress', hourglassControl='enhanced')
        
    def set_section_assignment(self):
        
        #* Stack direction of plate,
        #* the reference face is the top surface, the stacking direction is from bottom to top,
        #* The 1st ply is in the bottom surface (z0)
        #* The top face is the z1 face

        myPrt = self.model.parts[self.name_part]
        
        myPrt.SectionAssignment(region=myPrt.sets['all'], sectionName='orthotropic', offset=0.0, 
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
        myPrt.seedEdgeByNumber(edges=edges,
                    number=self.pMesh['num_element_thickness'], constraint=FIXED)

        #* Face edges
        self._seed_edge_face_hole_radial(myPrt, z0, reverse=False)
        self._seed_edge_face_circumferential_partition(myPrt, z0)
 
        if z1 == self.length_z:
            
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
        self.create_geometry_set('ply-%d'%(i_ply+1), points, geometry='cell')

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


class Bolt(Part):

    def __init__(self, model, name_part,
                length_z, r_hole, r_washer, thick_head,
                seedPart_size):
        
        super(Bolt,self).__init__(model, pGeo=None, pMesh=None)
        
        self.name_part = name_part
        self.length_z = length_z
        self.r_hole = r_hole
        self.r_washer = r_washer
        self.thick_head = thick_head
        self.seedPart_size = seedPart_size
        
        self.z0 = self.thick_head
        self.z1 = self.thick_head + self.length_z
        
    def create_sketch(self):
        '''
        Create sketch for bolts (side view).
        '''
        mySkt = self.model.ConstrainedSketch(name='bolt_sketch', sheetSize=200.0)
        
        mySkt.Line(point1=tuple([0.0, 0.0]),
                   point2=tuple([self.r_washer, 0.0]))
        mySkt.Line(point1=tuple([self.r_washer, 0.0]),
                   point2=tuple([self.r_washer, self.thick_head]))
        mySkt.Line(point1=tuple([self.r_washer, self.thick_head]),
                   point2=tuple([self.r_hole, self.thick_head]))
        mySkt.Line(point1=tuple([self.r_hole, self.thick_head]),
                   point2=tuple([self.r_hole, self.thick_head+self.length_z]))
        mySkt.Line(point1=tuple([self.r_hole, self.thick_head+self.length_z]),
                   point2=tuple([self.r_washer, self.thick_head+self.length_z]))
        mySkt.Line(point1=tuple([self.r_washer, self.thick_head+self.length_z]),
                   point2=tuple([self.r_washer, 2*self.thick_head+self.length_z]))
        mySkt.Line(point1=tuple([self.r_washer, 2*self.thick_head+self.length_z]),
                   point2=tuple([0.0, 2*self.thick_head+self.length_z]))
        mySkt.Line(point1=tuple([0.0, 2*self.thick_head+self.length_z]),
                   point2=tuple([0.0, 0.0]))

    def create_part(self):
        '''
        Use sketches for revolution
        '''
        myPrt = self.model.Part(name=self.name_part, dimensionality=THREE_D, type=DEFORMABLE_BODY)

        myPrt.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
        self.rename_feature(myPrt, 'XYPLANE')
        myPrt.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
        self.rename_feature(myPrt, 'XZPLANE')
        myPrt.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
        self.rename_feature(myPrt, 'YZPLANE')
        myPrt.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=self.z0)
        self.rename_feature(myPrt, 'XYPLANE-z0')
        myPrt.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=self.z1)
        self.rename_feature(myPrt, 'XYPLANE-z1')

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
            sketchPlane=self.get_datum_by_name(myPrt, 'YZPLANE'),
            sketchUpEdge=self.get_datum_by_name(myPrt, 'YAXIS'), 
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 0.0))
    
        #* Section sketch
        mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=200, transform=transform)
        mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)
        mySkt.retrieveSketch(sketch=self.model.sketches['bolt_sketch'])
        
        #* Part by revolution
        cl = mySkt.ConstructionLine(point1=(0.0, 0.0), point2=(0.0, 100.0))
        mySkt.assignCenterline(line=cl)
        
        myPrt.SolidRevolve(sketchPlane=self.get_datum_by_name(myPrt, 'YZPLANE'),
            sketchUpEdge=self.get_datum_by_name(myPrt, 'YAXIS'),
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM,
            sketch=mySkt, angle=360.0, flipRevolveDirection=OFF)

        #* Post procedure
        myPrt.setValues(geometryRefinement=EXTRA_FINE)
        del self.model.sketches['__profile__']

    def create_partition(self):

        myPrt = self.model.parts[self.name_part]
        myPrt.PartitionCellByDatumPlane(datumPlane=self.get_datum_by_name(myPrt, 'XYPLANE-z0'), cells=myPrt.cells)
        myPrt.PartitionCellByDatumPlane(datumPlane=self.get_datum_by_name(myPrt, 'XYPLANE-z1'), cells=myPrt.cells)
        myPrt.PartitionCellByDatumPlane(datumPlane=self.get_datum_by_name(myPrt, 'XZPLANE'), cells=myPrt.cells)
        myPrt.PartitionCellByDatumPlane(datumPlane=self.get_datum_by_name(myPrt, 'YZPLANE'), cells=myPrt.cells)
    
    def create_surface(self):
        myPrt = self.model.parts[self.name_part]

        #* Inner surfaces of bolt nut
        pts_z0 = []
        pts_z1 = []
        rr = 0.5*(self.r_hole+self.r_washer)
        for k in range(4):
            angle = (k+0.5)*np.pi/2.0
            pts_z0.append((rr*np.cos(angle), rr*np.sin(angle), self.z0))
            pts_z1.append((rr*np.cos(angle), rr*np.sin(angle), self.z1))
            
        faces = self.get_faces(myPrt, pts_z0, getClosest=True, searchTolerance=1e-3)
        myPrt.Surface(side1Faces=faces, name='bolt_head_inner_z0')
        faces = self.get_faces(myPrt, pts_z1, getClosest=True, searchTolerance=1e-3)
        myPrt.Surface(side1Faces=faces, name='bolt_head_inner_z1')
        
        #* Surface of bolt shank
        pts = []
        z = 0.5*(self.z0 + self.z1)
        for k in range(4):
            angle = (k+0.5)*np.pi/2.0
            pts.append((self.r_hole*np.cos(angle), self.r_hole*np.sin(angle), z))

        faces = self.get_faces(myPrt, pts, getClosest=True, searchTolerance=1e-3)
        myPrt.Surface(side1Faces=faces, name='bolt_shank')
    
    def create_set(self):
        myPrt = self.model.parts[self.name_part]
        myPrt.Set(cells=myPrt.cells, name='all')

    def set_seeding(self):
        myPrt = self.model.parts[self.name_part]
        myPrt.seedPart(size=self.seedPart_size, deviationFactor=0.1, minSizeFactor=0.1)

    def create_mesh(self):
        myPrt = self.model.parts[self.name_part]
        myPrt.generateMesh()
    
    def set_element_type(self):
        myPrt = self.model.parts[self.name_part]
        self.set_element_type_of_part(myPrt, kind='3D stress')

    def set_section_assignment(self):
        myPrt = self.model.parts[self.name_part]
        myPrt.SectionAssignment(region=myPrt.sets['all'], sectionName='Ti-6Al-4V', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


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
        
        self.create_material_titanium()
        self.create_section_titanium()

    def setup_parts(self):
        
        self.length_x = self.pGeo['len_x_plate']
        self.length_y = self.pGeo['len_y_plate']
        self.length_z = self.pGeo['len_z_plate']
        self.r_hole = self.pGeo['fasteners'][0]['r_hole']
        self.r_washer = self.pGeo['fasteners'][0]['r_washer']
        self.thick_head = self.pGeo['fasteners'][0]['thick_head']
        self.xc_hole = self.pGeo['fasteners'][0]['x_center']
        self.yc_hole = self.pGeo['fasteners'][0]['y_center']
        
        self.bolt = Bolt(self.model, name_part='bolt',
                    length_z=self.length_z*2,
                    r_hole=self.r_hole,
                    r_washer=self.r_washer,
                    thick_head=self.thick_head,
                    seedPart_size=self.pMesh['bolt_seedPart_size'])
        self.bolt.build()
        
        self.plate_0 = OpenHolePlate(self.model,
                            name_part='plate_0',
                            pMesh=self.pMesh,
                            length_x=self.length_x,
                            length_y=self.length_y,
                            length_z=self.length_z,
                            xc_hole=self.xc_hole,
                            yc_hole=self.yc_hole,
                            r_hole=self.r_hole)
        self.plate_0.build()
        
        self.plate_1 = OpenHolePlate(self.model,
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
        a.Instance(name='bolt', part=self.model.parts['bolt'], dependent=ON)
        a.Instance(name='plate_0', part=self.model.parts['plate_0'], dependent=ON)
        a.Instance(name='plate_1', part=self.model.parts['plate_1'], dependent=ON)

        a.translate(instanceList=('bolt',),
                    vector=(self.xc_hole, self.yc_hole, -self.thick_head))

        dx_plate1 = self.xc_hole - (self.length_x - self.xc_hole)
        a.translate(instanceList=('plate_1',),
                    vector=(dx_plate1, 0.0, self.length_z))

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
        '''
        Define interactions
        '''
        self.create_interaction_property_contact()
        
        self.model.ContactStd(name='General Interaction', createStepName='Initial')
        
        gi = self.model.interactions['General Interaction']
        gi.includedPairs.setValuesInStep(stepName='Initial', useAllstar=ON)
        gi.contactPropertyAssignments.appendInStep(stepName='Initial', 
                assignments=((GLOBAL, SELF, 'IP-Contact-Friction'), ))
        
    def setup_loads(self):

        a = self.rootAssembly
        
        self.create_reference_point(self.length_x*2, 0.5*self.length_y, 0, self.label_rp)
        self.create_reference_point_set(self.label_rp, self.label_rp)  

        #* Encastre BC on the left edge of plate_0
        self.model.EncastreBC(name='BC-x0', createStepName='Initial',
            region=a.sets['plate_0.face_x0'],
            localCsys=None)
        
        #* Coupling BC on the right edge of plate_1
        # Only constrain the loaded DOFs so the end face can deform freely in other directions
        _d = self.displacement
        self.model.Coupling(name='Coupling-plate_1_face_x1',
            controlPoint=a.sets[self.label_rp],
            surface=a.instances['plate_1'].surfaces['face_x1'],
            influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
            alpha=0.0, localCsys=None,
            u1=ON if _d[0] != 0 else OFF,
            u2=ON if _d[1] != 0 else OFF,
            u3=ON if _d[2] != 0 else OFF,
            ur1=OFF, ur2=OFF, ur3=OFF)
        
        # #* Apply displacement BCs on the reference points
        self.model.DisplacementBC(name=self.label_rp, createStepName='Loading', 
            region=a.sets[self.label_rp],
            u1=_d[0], u2=_d[1], u3=_d[2],
            ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF,
            distributionType=UNIFORM, fieldName='', localCsys=None)
        
    def setup_outputs(self):

        self.model.fieldOutputRequests['F-Output-1'].setValues(
            variables=('S', 'E', 'U', 'RF'), frequency=LAST_INCREMENT)
        
        self.model.FieldOutputRequest(name='FO-layup-0', 
            createStepName='Loading', variables=('S', 'E', 'SDV'), frequency=LAST_INCREMENT,
            layupNames=('plate_0.all', ), 
            layupLocationMethod=ALL_LOCATIONS, rebar=EXCLUDE)
        
        self.model.FieldOutputRequest(name='FO-layup-1', 
            createStepName='Loading', variables=('S', 'E', 'SDV'), frequency=LAST_INCREMENT,
            layupNames=('plate_1.all', ), 
            layupLocationMethod=ALL_LOCATIONS, rebar=EXCLUDE)
        
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


def extract_field(name_job, fname_save='specimen-field-C3D8R.dat'):
    
    N_SET = 2
    NAME_INSTANCES = ['PLATE_0', 'PLATE_1']
    NAME_SETS = ['PARTITION_CIRCLE', 'PARTITION_CIRCLE']
    
    odb = OdbOperation(name_job)

    with open(fname_save, 'w') as f:
        
        f.write('Variables= X Y Z index S11 S22 S33 S12 S13 S23\n')
        
        for i_set in range(N_SET):
            
            element_labels, indices_fieldOutput = odb.get_element_labels_and_indices(
                name_instance=NAME_INSTANCES[i_set], name_set=NAME_SETS[i_set])
            coordinates = odb.probe_element_center_coordinate(
                name_instance=NAME_INSTANCES[i_set], element_label=element_labels)
            values_S = odb.probe_element_values(variable='S', index_fieldOutput=indices_fieldOutput)
            n_element = len(indices_fieldOutput)
            
            f.write('zone T=" %s %s ", I= %d\n'%(NAME_INSTANCES[i_set], NAME_SETS[i_set], n_element))
            for i in range(n_element):
                for j in range(3):
                    f.write(' %14.6E'%(coordinates[i][j]))
                f.write(' %d'%(indices_fieldOutput[i]))
                for j in range(6):
                    f.write(' %14.6E'%(values_S[i][j]))
                f.write('\n')
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

    #* Build model

    name_job = 'Job_OHP_%d_%d'%(index_run, index_case)

    model = SingleLapBoltedJoint(name_job, pGeo, pMesh, pRun, displacement=displacement)
    model.build()
    model.set_view()
    model.save_cae('OHP_%d_%d.cae'%(index_run, index_case))
    
    if not parameters['not_run_job']:
        
        model.write_job_inp(model.name_job)
        model.submit_job(model.name_job, only_data_check=False)

        #* Post process
        odb = OdbOperation(model.name_job)

        with open(model.name_job+'-RF.dat', 'w') as f:

            rf_RP = odb.probe_node_values(variable='RF', index_fieldOutput=0)
            u_RP  = odb.probe_node_values(variable='U',  index_fieldOutput=0)

            for i in range(3):
                f.write('%s_RF%d  %20.6E \n'%(model.label_rp, i+1, rf_RP[i]))
            for i in range(3):
                f.write('%s_U%d   %20.6E \n'%(model.label_rp, i+1, u_RP[i]))

        extract_field(name_job=name_job, fname_save=name_job+'-field.dat')

