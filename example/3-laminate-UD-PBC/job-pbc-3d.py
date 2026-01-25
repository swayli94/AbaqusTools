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

import sys
import json

from laminate_C3D8R import LaminateModel
from AbaqusTools import OdbOperation
from AbaqusTools.pbc import PBC_3DOrthotropic


class TestModel_PBC_3D(LaminateModel):
    
    def __init__(self, name_job, pGeo, pMesh, pRun, strain_component=0, strain_scale=1E-3):
        
        super(TestModel_PBC_3D,self).__init__(name_job, pGeo, pMesh, pRun)

        self.strain_component = strain_component
        self.strain_scale = strain_scale
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
        name_instance = 'plate'
        
        #* Pairs of faces that are periodic
        #     master_face, slave_face,      forbidden_sets,     coords_sorting, name_mfn, name_sfn
        pairs = self._create_pbc_face_pairs()

        #* Create node sets on the master/slave faces
        label_forbidden = []

        for master_face, slave_face, forbidden_sets, coords_sorting, name_mfn, name_sfn in pairs:
            
            _,_, label_forbidden = PBC_3DOrthotropic.create_node_sets(
                                myMdl=self.model, name_instance=name_instance, 
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

    def _create_pbc_face_pairs(self):
        '''
        Create pairs of faces that are periodic.
        
        Notes
        ---------------
        Using forbidden sets to exclude duplicated constraints is not perfect.
        The correct way is to use the automatic approach of `label_forbidden_nodes` 
        in `exclude_forbidden_nodes_pbc` (in `pbc.py`).
        '''
        #     master_face, slave_face,      forbidden_sets,     coords_sorting, name_mfn, name_sfn
        
        if self.strain_component == 0 or self.strain_component == 3:
            # epsilon_11, gamma_23 (yz-face should be intact, i.e., not be forbidden)
            pairs=[
                ('face_x1', 'face_x0',  [], (1,2),  'MFn-X', 'SFn-X'),
                ('face_y1', 'face_y0',  [], (2,0),  'MFn-Y', 'SFn-Y'),
                ('face_z1', 'face_z0',  [], (0,1),  'MFn-Z', 'SFn-Z')]
                    
        elif self.strain_component == 1 or self.strain_component == 4:
            # epsilon_22, gamma_13 (zx-face should be intact, i.e., not be forbidden)
            pairs=[
                ('face_y1', 'face_y0',  [], (2,0),  'MFn-Y', 'SFn-Y'),
                ('face_z1', 'face_z0',  [], (0,1),  'MFn-Z', 'SFn-Z'),
                ('face_x1', 'face_x0',  [], (1,2),  'MFn-X', 'SFn-X')]
        
        else:
            # epsilon_33, gamma_12 (xy-face should be intact, i.e., not be forbidden)
            pairs=[
                ('face_z1', 'face_z0',  [], (0,1),  'MFn-Z', 'SFn-Z'),
                ('face_x1', 'face_x0',  [], (1,2),  'MFn-X', 'SFn-X'),
                ('face_y1', 'face_y0',  [], (2,0),  'MFn-Y', 'SFn-Y')]
            
        return pairs

    def create_bc_pinned(self):
        '''
        Create Abaqus BCs (After `Step`) to remove the rigid body motion.
        '''           
        a = self.rootAssembly
        
        self.model.DisplacementBC(name='Pinned-X0Y1Z0', createStepName='Initial', 
            region=a.sets['plate.vertex_010'],
            u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

        self.model.DisplacementBC(name='Pinned-X1Y1Z0', createStepName='Initial', 
            region=a.sets['plate.vertex_110'],
            u1=UNSET, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
        
        self.model.DisplacementBC(name='Pinned-X1Y0Z0', createStepName='Initial', 
            region=a.sets['plate.vertex_100'],
            u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
        
    def create_bc_displacement(self):
        '''
        Create Abaqus BCs (After `Step`) to apply the displacement boundary conditions.
        '''      
        a = self.rootAssembly

        for i_rp, label_rp in enumerate(self.label_rp):
            
            if i_rp == self.strain_component:
                u1 = self.strain_scale
            else:
                u1 = 0.0

            self.model.DisplacementBC(name=label_rp, createStepName='Loading', 
                region=a.sets[label_rp],
                u1=u1, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)


if __name__ == '__main__':

    with open('parameters.json', 'r') as f:
        parameters = json.load(f)

    pGeo = parameters['pGeo']
    pMesh = parameters['pMesh']
    pRun = parameters['pRun']
    
    index_run = parameters['index_run']
    index_strain_vector = parameters['index_strain_vector']

    print('>>> ')
    print('>>> Strain component: %d'%(index_strain_vector))
    print('>>> ')

    #* Build model

    name_job = 'Job_laminate_%d_%d'%(index_run, index_strain_vector)

    model = TestModel_PBC_3D(name_job, pGeo, pMesh, pRun,
                        strain_component=index_strain_vector,
                        strain_scale=parameters['strain_scale'])

    model.build()
    model.set_view()
    model.save_cae('laminate_%d_%d.cae'%(index_run, index_strain_vector))
    
    if not parameters['not_run_job']:
        
        model.write_job_inp(model.name_job)
        model.submit_job(model.name_job, only_data_check=False)
        
        #* Post process
        cmd_arguments = str(sys.argv)
        
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

                f.write('Strain_%d  %20.6E \n'%(index_strain_vector, model.strain_scale))
                