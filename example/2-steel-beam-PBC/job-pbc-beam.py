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

from steel_beam_C3D8R import SteelBeamModel
from AbaqusTools import OdbOperation
from AbaqusTools.pbc import PBC_Beam


class TestModel_PBC_Beam(SteelBeamModel):
    
    def __init__(self, name_job, strain_vector=[1E-6,0,0,0]):
        
        super(TestModel_PBC_Beam,self).__init__(name_job)

        self.strain_vector = strain_vector

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
        pairs = [ ('beam_0',           'face_z1',         'face_z0',                   [])]

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
            
            name_mfn, name_sfn, label_forbidden = PBC_Beam.create_node_sets(self.model, name_instance, 
                                        name_master_face, name_slave_face, coords_sorting, name_forbidden_sets, label_forbidden)
            
            self.face_pairs_name_node_sets.append((name_instance, name_mfn, name_sfn))

    def create_periodic_bc(self):
        '''
        Create Abaqus BCs for periodic boundary conditions
        '''
        #* Master node (reference point)
        self.create_reference_point(0, 0,  self.beam_0.length_z*1.1, 'myRP-1')
        self.create_reference_point_set('RP-1', 'myRP-1')
        
        self.create_reference_point(0, 20, self.beam_0.length_z*1.1, 'myRP-2')
        self.create_reference_point_set('RP-2', 'myRP-2')
        
        #* Create constraint equations
        for name_instance, name_mfn, name_sfn in self.face_pairs_name_node_sets:
            
            PBC_Beam.create_constraints_strain_vector(self.model, 'PBC_b-%s'%(name_instance), 
                            name_mfn, name_sfn, 'RP-1', 'RP-2', 
                            neutral_axis_x=self.neutral_axis_x, neutral_axis_y=self.neutral_axis_y)
    
        #* Define BCs (After `Step`)
        mp = self.beam_0
        
        a = self.rootAssembly
        
        self.model.DisplacementBC(name='Pinned-X0Y1Z0', createStepName='Initial', 
            region=a.sets['beam_0.vertex_010'],
            u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

        self.model.DisplacementBC(name='Pinned-X1Y1Z0', createStepName='Initial', 
            region=a.sets['beam_0.vertex_110'],
            u1=UNSET, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
        
        self.model.DisplacementBC(name='Pinned-X1Y0Z0', createStepName='Initial', 
            region=a.sets['beam_0.vertex_100'],
            u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
        
        
        u3_MN1, u1_MN2, u2_MN2, u3_MN2 = PBC_Beam.calculate_master_node_displacement_BC(self.strain_vector, self.beam_0.length_z)
        
        self.model.DisplacementBC(name='RP-1', createStepName='Loading', 
            region=a.sets['RP-1'],
            u1=UNSET, u2=UNSET, u3=u3_MN1, 
            ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
            fieldName='', localCsys=None)
        
        self.model.DisplacementBC(name='RP-2', createStepName='Loading', 
            region=a.sets['RP-2'],
            u1=u1_MN2, u2=u2_MN2, u3=u3_MN2, 
            ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
            fieldName='', localCsys=None)
        

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

    model = TestModel_PBC_Beam(name_job, STRAIN_VECTORS[i0,:])
    model.build()
    
    model.write_job_inp(model.name_job)
    model.save_cae('beam_%d.cae'%(i0))
    
    model.submit_job(model.name_job, only_data_check=False)
    model.set_view()

    #* Post process
    if '-noGUI' in cmd_arguments:
    
        odb = OdbOperation(model.name_job)
        rf_RP1 = odb.probe_node_values(variable='RF', index_fieldOutput=0)
        rf_RP2 = odb.probe_node_values(variable='RF', index_fieldOutput=1)
        u_RP1  = odb.probe_node_values(variable='U',  index_fieldOutput=0)
        u_RP2  = odb.probe_node_values(variable='U',  index_fieldOutput=1)

        with open(model.name_job+'-RF.dat', 'w') as f:
            
            f.write('RP-1_RF3     %20.6E \n'%(rf_RP1[2]))
            f.write('RP-2_RF1     %20.6E \n'%(rf_RP2[0]))
            f.write('RP-2_RF2     %20.6E \n'%(rf_RP2[1]))
            f.write('RP-2_RF3     %20.6E \n'%(rf_RP2[2]))
            f.write('\n')
            f.write('RP-1_U3      %20.6E \n'%(u_RP1[2]))
            f.write('RP-2_U1      %20.6E \n'%(u_RP2[0]))
            f.write('RP-2_U2      %20.6E \n'%(u_RP2[1]))
            f.write('RP-2_U3      %20.6E \n'%(u_RP2[2]))
            f.write('\n')
            f.write('Strain_%d-1  %20.6E \n'%(i0, STRAIN_VECTORS[i0][0]))
            f.write('Strain_%d-2  %20.6E \n'%(i0, STRAIN_VECTORS[i0][1]))
            f.write('Strain_%d-3  %20.6E \n'%(i0, STRAIN_VECTORS[i0][2]))
            f.write('Strain_%d-4  %20.6E \n'%(i0, STRAIN_VECTORS[i0][3]))
            f.write('\n')
            f.write('RP-1_U1      %20.6E \n'%(u_RP1[0]))
            f.write('RP-1_U2      %20.6E \n'%(u_RP1[1]))
            f.write('\n')
            
