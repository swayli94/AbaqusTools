'''
A laminate with 3D stress elements (C3D8R).

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
from AbaqusTools.functions import LayupParameters


class Plate(Part):
    '''
    Plate with an open hole.
    '''
    def __init__(self, model, pGeo, pMesh):
        
        super(Plate, self).__init__(model, pGeo, pMesh)
        
        self.name_part = 'plate'

        #* Attributes
        self.length_x = self.pGeo['len_x_plate']
        self.length_y = self.pGeo['len_y_plate']
        self.length_z = self.pGeo['len_z_plate']
        
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
            mySkt.Line(point1=tuple(sketch_points[i-1,:]), point2=tuple(sketch_points[i,:]))

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
        myPrt.BaseSolidExtrude(sketch=mySkt, depth=self.length_z)

        #* Post procedure
        myPrt.setValues(geometryRefinement=EXTRA_FINE)
        del self.model.sketches['__profile__']
    
    #* Surface, set for the entire part
    
    def create_partition(self):
        pass
     
    def create_surface(self):
        
        myPrt = self.model.parts[self.name_part]

        faces = self.get_faces(myPrt, (0, 0.5*self.length_y, 0.5*self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_x0')
        
        faces = self.get_faces(myPrt, (self.length_x, 0.5*self.length_y, 0.5*self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_x1')

        faces = self.get_faces(myPrt, (0.5*self.length_x, 0.0, 0.5*self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_y0')
        
        faces = self.get_faces(myPrt, (0.5*self.length_x, self.length_y, 0.5*self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_y1')
        
        faces = self.get_faces(myPrt, (0.5*self.length_x, 0.5*self.length_y, 0.0))
        myPrt.Surface(side1Faces=faces, name='face_z0')
        
        faces = self.get_faces(myPrt, (0.5*self.length_x, 0.5*self.length_y, self.length_z))
        myPrt.Surface(side1Faces=faces, name='face_z1')

    def create_set(self):
        
        myPrt = self.model.parts[self.name_part]
        myPrt.Set(cells=myPrt.cells, name='all') 

        lx = self.length_x
        ly = self.length_y
        lz = self.length_z
        
        self.create_geometry_set('face_x0', (0.0,    0.5*ly, 0.5*lz), geometry='face')
        self.create_geometry_set('face_x1', (lx,     0.5*ly, 0.5*lz), geometry='face')
        self.create_geometry_set('face_y0', (0.5*lx, 0.0,    0.5*lz), geometry='face')
        self.create_geometry_set('face_y1', (0.5*lx, ly,     0.5*lz), geometry='face')
        self.create_geometry_set('face_z0', (0.5*lx, 0.5*ly, 0.0   ), geometry='face')
        self.create_geometry_set('face_z1', (0.5*lx, 0.5*ly, lz    ), geometry='face')

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
        
        self.create_geometry_set('vertex_000', (0.0, 0.0, 0.0), geometry='vertex')
        self.create_geometry_set('vertex_100', (lx,  0.0, 0.0), geometry='vertex')
        self.create_geometry_set('vertex_010', (0.0, ly,  0.0), geometry='vertex')
        self.create_geometry_set('vertex_110', (lx,  ly,  0.0), geometry='vertex')
        self.create_geometry_set('vertex_001', (0.0, 0.0, lz ), geometry='vertex')
        self.create_geometry_set('vertex_101', (lx,  0.0, lz ), geometry='vertex')
        self.create_geometry_set('vertex_011', (0.0, ly,  lz ), geometry='vertex')
        self.create_geometry_set('vertex_111', (lx,  ly,  lz ), geometry='vertex')
    

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
        points=[(0.5*self.length_x, 0.5*self.length_y, z_mid)]
        self.create_geometry_set('ply-%d'%(i_ply+1), points, geometry='cell')


class LaminateModel(Model):
    
    def __init__(self, name_job, pGeo, pMesh, pRun):
        
        super(LaminateModel,self).__init__(pGeo=pGeo, pMesh=pMesh, pRun=pRun, name_job=name_job)
        
    def initialization(self):
        
        self.model = mdb.models[str(self.name_model)]

        self.create_material_IM785517(elastic_type='ENGINEERING_CONSTANTS')
        
        self.create_section_IM785517()

    def setup_parts(self):
        
        self.plate = Plate(self.model, self.pGeo, self.pMesh)
        self.plate.build()
        
        self.length_x = self.plate.length_x
        self.length_y = self.plate.length_y
        self.length_z = self.plate.length_z
        self.volume = self.length_x*self.length_y*self.length_z
    
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

        self.model.fieldOutputRequests['F-Output-1'].setValues(
            variables=('S', 'E', 'U', 'RF'), frequency=LAST_INCREMENT)
        
        self.model.FieldOutputRequest(name='FO-layup', 
            createStepName='Loading', variables=('S', 'E', 'SDV'), frequency=LAST_INCREMENT,
            layupNames=('plate.all', ), 
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


if __name__ == '__main__':

    with open('default-parameters.json', 'r') as f:
        default_parameters = json.load(f)

    pGeo = default_parameters['pGeo']
    pMesh = default_parameters['pMesh']
    pRun = default_parameters['pRun']

    model = LaminateModel('Job_LMT', pGeo, pMesh, pRun)
    
    model.build()
    model.set_view()
