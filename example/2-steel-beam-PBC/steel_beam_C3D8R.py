'''
A simple steel beam with 3D stress elements (C3D8R).


Example
---------------
Put the python script in the same directory as the 'AbaqusTools' folder.
Clean up the `*.pyc` files first, which are compiled by Abaqus in previous jobs.
The `clean.bat` gives an example of cleaning.

>>> abaqus cae script=*.py
'''
from abaqus import *
from abaqusConstants import *
from sketch import *

from AbaqusTools import Part, Model


class Beam(Part):

    def __init__(self, model, name_part, length_x=100, length_y=20, length_z=150, 
                    seedPart_size=5):
        
        super(Beam,self).__init__(model, None, None)
        
        self.name_part = name_part
        
        self.length_x = length_x
        self.length_y = length_y
        self.length_z = length_z
        
        self.seedPart_size = seedPart_size

    def create_sketch(self):

        mySkt = self.model.ConstrainedSketch(name='side_view', sheetSize=200)
        mySkt.Line(point1=(0, 0),               point2=(self.length_x, 0))
        mySkt.Line(point1=(0, 0),               point2=(0,             self.length_y))
        mySkt.Line(point1=(0, self.length_y),   point2=(self.length_x, self.length_y))
        mySkt.Line(point1=(self.length_x, 0),   point2=(self.length_x, self.length_y))

    def create_part(self):

        myPrt = self.model.Part(name=self.name_part, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        
        #* Reference plane and axis
        myPrt.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
        myPrt.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
        myPrt.DatumAxisByPrincipalAxis(principalAxis=XAXIS)
        datums = myPrt.datums
        
        #* Plane for wing sketch
        t1 = myPrt.MakeSketchTransform(sketchPlane=datums[1], sketchUpEdge=datums[3], 
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 0.0))
        
        #* Section sketch
        mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=200, transform=t1)
        mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)
        mySkt.retrieveSketch(sketch=self.model.sketches['side_view'])

        #* Part by extrusion
        myPrt.BaseSolidExtrude(sketch=mySkt, depth=self.length_z)
        del self.model.sketches['__profile__']

        #* Post procedure
        myPrt.setValues(geometryRefinement=EXTRA_FINE)
        
    def create_surface(self):
    
        myPrt = self.model.parts[self.name_part]

        faces = self.get_faces(myPrt, (0.5*self.length_x, self.length_y, 0.01*self.length_z))
        myPrt.Surface(side1Faces=faces, name='y1')
        
        faces = self.get_faces(myPrt, (0.5*self.length_x, 0.5*self.length_y, self.length_z))
        myPrt.Surface(side1Faces=faces, name='z1')

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

    def set_seeding(self):

        myPrt = self.model.parts[self.name_part]
        myPrt.seedPart(size=self.seedPart_size, deviationFactor=0.1, minSizeFactor=0.1)
        
    def create_mesh(self):

        myPrt = self.model.parts[self.name_part]
        myPrt.setMeshControls(regions=myPrt.cells, technique=STRUCTURED)
        myPrt.assignStackDirection(referenceRegion=myPrt.surfaces['y1'].faces[0], cells=myPrt.cells)
        myPrt.generateMesh()
    
    def set_element_type(self):

        myPrt = self.model.parts[self.name_part]
        self.set_element_type_of_part(myPrt, kind='3D stress')

    def set_section_assignment(self):
        
        myPrt = self.model.parts[self.name_part]
        
        myPrt.SectionAssignment(region=myPrt.sets['all'], sectionName='Steel', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


class SteelBeamModel(Model):
    
    def __init__(self, name_job):
        
        super(SteelBeamModel,self).__init__(pGeo=None, pMesh=None, pRun=None)

        self.name_job = name_job
        
    def initialization(self):
        
        self.model = mdb.models[str(self.name_model)]

        self.create_material_steel()

        self.create_section_steel()

    def setup_parts(self):
        
        self.beam_0 = Beam(self.model, 'beam_0', 
                        length_x=10, length_y=10, length_z=10,
                        seedPart_size=1.0)
        self.beam_0.build()
        
        self.neutral_axis_x = 0.5*self.beam_0.length_x
        self.neutral_axis_y = 0.5*self.beam_0.length_y
        
        self.length_x = self.beam_0.length_x
        self.length_y = self.beam_0.length_y
        self.length_z = self.beam_0.length_z
        self.volume = self.length_x*self.length_y*self.length_z
    
    def setup_assembly(self):
        
        a = self.rootAssembly
        p = self.model.parts['beam_0']
        a.Instance(name='beam_0', part=p, dependent=ON)
    
    def setup_steps(self):

        self.create_static_step(timePeriod=1.0, maxNumInc=10000, 
                initialInc=0.01, minInc=1E-15, maxInc=0.1, nlgeom=False)
        
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

        mdb.Job(name=self.name_job, model=str(self.name_model), description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
            multiprocessingMode=DEFAULT, numCpus=4, numDomains=4, numGPUs=0)


if __name__ == '__main__':

    name_job = 'Demo_steel_beam_C3D8R'

    model = SteelBeamModel(name_job)
    model.build()
    model.set_view()
