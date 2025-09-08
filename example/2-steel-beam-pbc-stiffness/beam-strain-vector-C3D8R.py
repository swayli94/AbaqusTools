'''
Test the Periodic Boundary Condition (PBC) in Abaqus with a simple steel beam.

Define PBC by constraint equations between node sets.

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

import os
import sys
import numpy as np

from AbaqusTools import Part, Model, OdbOperation
from AbaqusTools.pbc import PeriodicBC


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

        faces = self.get_faces(myPrt, (0.5*self.length_x, 0.5*self.length_y, 0))
        myPrt.Set(faces=faces, name='z0')

        faces = self.get_faces(myPrt, (0.5*self.length_x, 0.5*self.length_y, self.length_z))
        myPrt.Set(faces=faces, name='z1')

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
    
    def __init__(self, name_job, strain_vector=[1E-6,0,0,0]):
        
        super(TestModel,self).__init__(pGeo=None, pMesh=None, pRun=None)

        self.name_job = name_job
        
        self.strain_vector = strain_vector

    def initialization(self):
        
        self.model = mdb.models[str(self.name_model)]

        self.create_material_steel()

        self.create_section_steel()

    def setup_parts(self):
        
        self.beam_0 = Beam(self.model, 'beam_0', length_x=100, length_y=20, length_z=50, seedPart_size=4.0)
        self.beam_0.build()
        
        self.neutral_axis_x = 0.5*self.beam_0.length_x
        self.neutral_axis_y = 0.5*self.beam_0.length_y
    
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
        self.setup_periodic_bc_node_sets()

        self.create_periodic_bc()
    
    def setup_periodic_bc_node_sets(self):
        '''
        Create node sets of nodes in paired faces.
        '''
        #* Pairs of faces that are periodic
        #        name_instance,     name_master_face,   name_slave_face,    name_forbidden_sets
        pairs = [ ('beam_0',           'z1',               'z0',                   [])]

        # (name_instance, name_master_face_node_set, name_slave_face_node_set)
        self.face_pairs_name_node_sets = []

        #* Create node sets on the master/slave faces
        coords_sorting = (0, 1)
        last_instance = None
        
        for name_instance, name_master_face, name_slave_face, name_forbidden_sets in pairs:

            if last_instance is None:
                label_forbidden = []
            elif last_instance != name_instance:
                label_forbidden = []
            last_instance = name_instance
            
            name_mfn, name_sfn = PeriodicBC.create_node_sets(self.model, name_instance, 
                                        name_master_face, name_slave_face, coords_sorting, name_forbidden_sets, label_forbidden)
            
            self.face_pairs_name_node_sets.append((name_instance, name_mfn, name_sfn))

    def create_periodic_bc(self):
        '''
        Create Abaqus BCs for periodic boundary conditions
        '''
        #* Master node (reference point)
        self.create_reference_point(0, 0,  self.beam_0.length_z*1.1, 'MasterNode-1')
        self.create_reference_point_set('MasterNode-1', 'MasterNode-1')
        
        self.create_reference_point(0, 20, self.beam_0.length_z*1.1, 'MasterNode-2')
        self.create_reference_point_set('MasterNode-2', 'MasterNode-2')
        
        #* Create constraint equations
        for name_instance, name_mfn, name_sfn in self.face_pairs_name_node_sets:
            
            PeriodicBC.create_constraints_strain_vector(self.model, 'PBC_b-%s'%(name_instance), 
                            name_mfn, name_sfn, 'MasterNode-1', 'MasterNode-2', 
                            neutral_axis_x=self.neutral_axis_x, neutral_axis_y=self.neutral_axis_y)
    
        #* Define BCs (After `Step`)
        mp = self.beam_0
        
        pt = (0.0, mp.length_y, 0)
        mp.create_geometry_set('X0Y1Z0', 'beam_0', pt, geometry='vertex', getClosest=False)
        
        pt = (mp.length_x, 0.0, 0)
        mp.create_geometry_set('X1Y0Z0', 'beam_0', pt, geometry='vertex', getClosest=False)
        
        pt = (mp.length_x, mp.length_y, 0)
        mp.create_geometry_set('X1Y1Z0', 'beam_0', pt, geometry='vertex', getClosest=False)
        
        a = self.rootAssembly
        
        self.model.DisplacementBC(name='Pinned-X0Y1Z0', createStepName='Initial', 
            region=a.sets['beam_0.X0Y1Z0'],
            u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

        self.model.DisplacementBC(name='Pinned-X1Y1Z0', createStepName='Initial', 
            region=a.sets['beam_0.X1Y1Z0'],
            u1=UNSET, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
        
        self.model.DisplacementBC(name='Pinned-X1Y0Z0', createStepName='Initial', 
            region=a.sets['beam_0.X1Y0Z0'],
            u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
        
        
        u3_MN1, u1_MN2, u2_MN2, u3_MN2 = PeriodicBC.calculate_master_node_displacement_BC(self.strain_vector, self.beam_0.length_z)
        
        self.model.DisplacementBC(name='MasterNode1', createStepName='Loading', 
            region=a.sets['MasterNode-1'],
            u1=UNSET, u2=UNSET, u3=u3_MN1, 
            ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
            fieldName='', localCsys=None)
        
        self.model.DisplacementBC(name='MasterNode2', createStepName='Loading', 
            region=a.sets['MasterNode-2'],
            u1=u1_MN2, u2=u2_MN2, u3=u3_MN2, 
            ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
            fieldName='', localCsys=None)
        
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


    cmd_arguments = str(sys.argv)

    STRAIN_VECTORS=np.eye(4,4)*1E-6

    #* Read strain vector from file

    if os.path.exists('temp-strain-vector.txt'):
        
        with open('temp-strain-vector.txt', 'r') as f:
            lines = f.readlines()
            i0 = int(lines[0].split()[0])

    else:
        
        i0 = 0
            
    print('>>> ')
    print('>>> Base strain vector: %d'%(i0))
    print('>>> ')

    name_job = 'Job_beam_%d'%(i0)

    model = TestModel(name_job, STRAIN_VECTORS[i0,:])
    model.build()
    
    model.write_job_inp(model.name_job)
    model.save_cae('beam_%d.cae'%(i0))
    
    model.submit_job(model.name_job, only_data_check=False)
    model.set_view()

    #* Post process
    if '-noGUI' in cmd_arguments:
    
        odb = OdbOperation(model.name_job)
        rf_MN1 = odb.probe_node_values(variable='RF', index_fieldOutput=0)
        rf_MN2 = odb.probe_node_values(variable='RF', index_fieldOutput=1)
        u_MN1  = odb.probe_node_values(variable='U',  index_fieldOutput=0)
        u_MN2  = odb.probe_node_values(variable='U',  index_fieldOutput=1)
        u_MN5  = odb.probe_node_values(variable='U',  index_fieldOutput=2)


        with open(model.name_job+'-RF.dat', 'w') as f:
            
            f.write('MN1_RF3     %20.6E \n'%(rf_MN1[2]))
            f.write('MN2_RF1     %20.6E \n'%(rf_MN2[0]))
            f.write('MN2_RF2     %20.6E \n'%(rf_MN2[1]))
            f.write('MN2_RF3     %20.6E \n'%(rf_MN2[2]))
            f.write('\n')
            f.write('MN1_U3      %20.6E \n'%(u_MN1[2]))
            f.write('MN2_U1      %20.6E \n'%(u_MN2[0]))
            f.write('MN2_U2      %20.6E \n'%(u_MN2[1]))
            f.write('MN2_U3      %20.6E \n'%(u_MN2[2]))
            f.write('\n')
            f.write('Strain_%d-1  %20.6E \n'%(i0, STRAIN_VECTORS[i0][0]))
            f.write('Strain_%d-2  %20.6E \n'%(i0, STRAIN_VECTORS[i0][1]))
            f.write('Strain_%d-3  %20.6E \n'%(i0, STRAIN_VECTORS[i0][2]))
            f.write('Strain_%d-4  %20.6E \n'%(i0, STRAIN_VECTORS[i0][3]))
            f.write('\n')
            f.write('MN1_U1      %20.6E \n'%(u_MN1[0]))
            f.write('MN1_U2      %20.6E \n'%(u_MN1[1]))
            f.write('\n')
            
