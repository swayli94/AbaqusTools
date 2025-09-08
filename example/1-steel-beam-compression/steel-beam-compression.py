'''
Test the `Part` and `Model` classes with a simple steel beam compression case.

Example
---------------
Put the python script in the same directory as the 'AbaqusTools' folder.
Clean up the `*.pyc` files first, which are compiled by Abaqus in previous jobs.
The `clean.bat` gives an example of cleaning.

>>> abaqus cae script=steel-beam-compression.py
>>> abaqus cae noGUI=steel-beam-compression.py
'''
from abaqus import *
from abaqusConstants import *
from sketch import *
import mesh

from AbaqusTools import Part, Model


class SteelBeam(Part):

    def __init__(self, model, name_part, length_x=100, length_y=20, length_z=150, 
                    seedPart_size=5):
        
        super(SteelBeam,self).__init__(model, None, None)
        
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

        faces = self.get_faces(myPrt, (0.5*self.length_x, 0.5*self.length_y, 0))
        myPrt.Set(faces=faces, name='z0')

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


class TestModel(Model):
    
    def __init__(self):
        
        super(TestModel,self).__init__(pGeo=None, pMesh=None, pRun=None)

    def initialization(self):
        
        self.model = mdb.models[str(self.name_model)]

        self.create_material_steel()

        self.create_section_steel()

    def setup_parts(self):
        
        self.beam_0 = SteelBeam(self.model, 'beam_0', length_x=100, length_y=20, length_z=150)
        self.beam_0.build()
    
    def setup_assembly(self):
        
        a = self.rootAssembly
        p = self.model.parts['beam_0']
        a.Instance(name='beam_0', part=p, dependent=ON)
    
    def setup_steps(self):

        self.create_static_step(timePeriod=1.0, maxNumInc=10000, 
                initialInc=0.01, minInc=1E-15, maxInc=0.1, nlgeom=False)
        
    def setup_loads(self):
        
        region = self.instance('beam_0').sets['z0']
        self.model.DisplacementBC(name='BC-wall-z0', createStepName='Initial', 
            region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
        
        region = self.instance('beam_0').surfaces['z1']
        self.model.Pressure(name='Load-pressure-z1', createStepName='Loading', 
            region=region, distributionType=UNIFORM, field='', magnitude=10, 
            amplitude=UNSET)
        
    def setup_outputs(self):

        self.model.fieldOutputRequests['F-Output-1'].setValues(
            variables=('S', 'E', 'U', 'RF'), 
            frequency=1)
        
    def setup_jobs(self):

        mdb.Job(name='Job-1', model=str(self.name_model), description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
            multiprocessingMode=DEFAULT, numCpus=4, numDomains=4, numGPUs=0)


if __name__ == '__main__':

    model = TestModel()
    
    model.build()

    model.submit_job('Job-1', only_data_check=False)
    
    model.set_view()
