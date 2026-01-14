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
from AbaqusTools.pbc import PBC_3DOrthotropic


class TestModel_PBC_3D(SteelBeamModel):
    
    def __init__(self, name_job, strain_vector=[1E-6,0,0,0,0,0]):
        
        super(TestModel_PBC_3D,self).__init__(name_job)

        self.strain_vector = strain_vector
        self.label_rp = ['RP_11', 'RP_22', 'RP_33', 'RP_23', 'RP_13', 'RP_12']

    def setup_loads(self):
        '''
        Define loads and boundary conditions
        '''
        self.create_reference_points()

        self.create_pbc_constraints()
        
        # self.create_bc_pinned()
        
        self.create_bc_displacement()
    
    def create_reference_points(self):
        '''
        Create reference points for the periodic boundary conditions
        '''
        for i_rp, label_rp in enumerate(self.label_rp):
            self.create_reference_point(-10*(i_rp+1), 0, 0, label_rp)
            self.create_reference_point_set(label_rp, label_rp)        
    
    def create_pbc_constraints(self):
        '''
        Create node sets of nodes in paired faces.
        '''
        name_instance = 'beam_0'
        
        #* Pairs of faces that are periodic
        #     master_face, slave_face,      forbidden_sets,     coords_sorting, name_mfn, name_sfn
        pairs=[('face_z1',  'face_z0',  [],                              (0,1),  'MFn-Z', 'SFn-Z'),
               ('face_x1',  'face_x0',  ['edge_y_z0x0', 'edge_y_z1x0', 
                                        'edge_y_z0x1', 'edge_y_z1x1'],   (1,2),  'MFn-X', 'SFn-X'),
               ('face_y1',  'face_y0',  ['edge_x_y0z0', 'edge_x_y1z0',
                                        'edge_x_y0z1', 'edge_x_y1z1',
                                        'edge_z_x0y0', 'edge_z_x1y0',
                                        'edge_z_x0y1', 'edge_z_x1y1'],   (2,0),  'MFn-Y', 'SFn-Y')]

        #* Create node sets on the master/slave faces
        label_forbidden = []

        for master_face, slave_face, forbidden_sets, coords_sorting, name_mfn, name_sfn in pairs:
            
            PBC_3DOrthotropic.create_node_sets(self.model, name_instance, 
                                name_master_face_set=master_face, 
                                name_slave_face_set=slave_face,
                                coords_sorting=coords_sorting,
                                name_forbidden_sets=forbidden_sets, 
                                label_forbidden_nodes=label_forbidden,
                                name_mfn=name_mfn, name_sfn=name_sfn)
        
        #* Create constraint equations
        PBC_3DOrthotropic.create_constraints_strain_vector(self.model,
                name_eqn='PBC3D_%s'%(name_instance),
                name_mfn_x_set='MFn-X', name_sfn_x_set='SFn-X',
                name_mfn_y_set='MFn-Y', name_sfn_y_set='SFn-Y',
                name_mfn_z_set='MFn-Z', name_sfn_z_set='SFn-Z',
                length_x=self.length_x, length_y=self.length_y, length_z=self.length_z,
                name_rp11=self.label_rp[0], name_rp22=self.label_rp[1], name_rp33=self.label_rp[2],
                name_rp23=self.label_rp[3], name_rp13=self.label_rp[4], name_rp12=self.label_rp[5])

    def create_bc_pinned(self):
        '''
        Create Abaqus BCs (After `Step`) to remove the rigid body motion.
        '''           
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
        
    def create_bc_displacement(self):
        '''
        Create Abaqus BCs (After `Step`) to apply the displacement boundary conditions.
        '''      
        a = self.rootAssembly

        for i_rp, label_rp in enumerate(self.label_rp):
        
            self.model.DisplacementBC(name=label_rp, createStepName='Loading', 
                region=a.sets[label_rp],
                u1=self.strain_vector[i_rp], 
                u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)


if __name__ == '__main__':


    cmd_arguments = str(sys.argv)

    STRAIN_VECTORS=np.eye(6,6)*1E-6

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

    model = TestModel_PBC_3D(name_job, STRAIN_VECTORS[i0,:])
    model.build()
    
    model.write_job_inp(model.name_job)
    model.save_cae('beam_%d.cae'%(i0))
    
    model.submit_job(model.name_job, only_data_check=False)
    model.set_view()

    #* Post process
    if '-noGUI' in cmd_arguments:
    
        odb = OdbOperation(model.name_job)

        with open(model.name_job+'-RF.dat', 'w') as f:
            
            rf_RPs = []
            u_RPs  = []
            
            for i_rp, label_rp in enumerate(model.label_rp):
                
                rf_RP = odb.probe_node_values(variable='RF', index_fieldOutput=i_rp)
                u_RP  = odb.probe_node_values(variable='U',  index_fieldOutput=i_rp)

                rf_RPs.append(rf_RP[0])
                u_RPs.append(u_RP[0])
            
            '''
            Note:
            
            The displacement `u1` of reference points is the strain `epsilon_ij`.
            Therefore, the stiffness components are the reaction forces divided by
            the volume of the model, instead of the area of the corresponding face.
            '''
                 
            for i_rp, label_rp in enumerate(model.label_rp):
                f.write('%s_RF  %20.6E \n'%(label_rp, rf_RPs[i_rp]/model.volume))
                
            for i_rp, label_rp in enumerate(model.label_rp):
                f.write('%s_U   %20.6E \n'%(label_rp, u_RPs[i_rp]))

            for i_rp, label_rp in enumerate(model.label_rp):
                f.write('Strain_%d-%d  %20.6E \n'%(i0, i_rp, STRAIN_VECTORS[i0][i_rp]))
            