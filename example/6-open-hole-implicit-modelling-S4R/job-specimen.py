'''
RVE with Periodic Boundary Condition (PBC).
(8_common)
'''
from abaqus import *
from abaqusConstants import *
from sketch import *

import sys
import json

from laminate_S4R import LaminateModel
from AbaqusTools import OdbOperation
from AbaqusTools.lin_bc import LBC_InPlaneLoad


class TestModel_InPlaneLoad(LaminateModel):
    
    def __init__(self, name_job, pGeo, pMesh, pRun, strain_component=0, strain_scale=1E-3):
        
        super(TestModel_InPlaneLoad,self).__init__(name_job, pGeo, pMesh, pRun)

        self.strain_component = strain_component
        self.strain_scale = strain_scale
        self.label_rp = ['RP_11', 'RP_22', 'RP_12']

    def setup_loads(self):
        '''
        Define loads and boundary conditions
        '''
        self.create_reference_points()

        self.create_lbc_constraints()

        self.create_bc_displacement()
    
    def create_reference_points(self):
        '''
        Create reference points for the boundary conditions
        '''
        for i_rp, label_rp in enumerate(self.label_rp):
            self.create_reference_point(-10*(i_rp+1), 0, 0, label_rp)
            self.create_reference_point_set(label_rp, label_rp)        

    def create_lbc_constraints(self):
        '''
        Create in-plane load boundary conditions using constraint equations.
        '''
        name_instance = 'plate'
        
        #* Pairs of edges
        if self.strain_component == 0:
            pairs=[ ('edge_x1', 'edge_x0',  [], (1,2),  'MFn-X', 'SFn-X')]
        elif self.strain_component == 1:
            pairs=[ ('edge_y1', 'edge_y0',  [], (2,0),  'MFn-Y', 'SFn-Y')]
        elif self.strain_component == 2:
            pairs=[ ('edge_x1', 'edge_x0',  [], (1,2),  'MFn-X', 'SFn-X'),
                    ('edge_y1', 'edge_y0',  [], (2,0),  'MFn-Y', 'SFn-Y')]
        else:
            raise ValueError('Invalid strain component index: %d'%(self.strain_component))

        #* Create node sets on the master/slave edges
        label_forbidden = []

        for master_face, slave_face, forbidden_sets, coords_sorting, name_mfn, name_sfn in pairs:
            
            _,_, label_forbidden = LBC_InPlaneLoad.create_node_sets(
                                myMdl=self.model, name_instance=name_instance, 
                                name_master_face_set=master_face, 
                                name_slave_face_set=slave_face,
                                coords_sorting=coords_sorting,
                                name_forbidden_sets=forbidden_sets, 
                                label_forbidden_nodes=label_forbidden,
                                name_mfn=name_mfn, name_sfn=name_sfn)
        
        #* Create constraint equations
        if self.strain_component == 0:
            LBC_InPlaneLoad.create_constraints_epsilon11(self.model,
                    name_eqn='InPlane11_%s'%(name_instance),
                    name_mfn_x_set='MFn-X',
                    length_x=self.length_x,
                    name_rp11=self.label_rp[0])
            
        elif self.strain_component == 1:
            LBC_InPlaneLoad.create_constraints_epsilon22(self.model,
                    name_eqn='InPlane22_%s'%(name_instance),
                    name_mfn_y_set='MFn-Y',
                    length_y=self.length_y, 
                    name_rp22=self.label_rp[1])
            
        elif self.strain_component == 2:
            LBC_InPlaneLoad.create_constraints_gamma12(self.model,
                    name_eqn='InPlane12_%s'%(name_instance),
                    name_mfn_x_set='MFn-X', name_sfn_x_set='SFn-X',
                    name_mfn_y_set='MFn-Y', name_sfn_y_set='SFn-Y',
                    length_x=self.length_x, length_y=self.length_y, 
                    name_rp12=self.label_rp[2])
        
    def create_bc_displacement(self):
        '''
        Create Abaqus BCs (After `Step`) to apply the displacement boundary conditions.
        '''      
        a = self.rootAssembly
        label_rp = self.label_rp[self.strain_component]

        #* Apply pinning BCs on the faces/edges/vertices to prevent rigid body motions
        if self.strain_component == 0:
            # epsilon_11
            self.model.EncastreBC(name='BC-x0', createStepName='Initial',
                region=a.sets['plate.edge_x0'],
                localCsys=None)
            
        elif self.strain_component == 1:
            # epsilon_22
            self.model.EncastreBC(name='BC-y0', createStepName='Initial',
                region=a.sets['plate.edge_y0'],
                localCsys=None)

        elif self.strain_component == 2:
            # gamma_12
            self.model.DisplacementBC(name='Pinned-X0Y0', createStepName='Initial', 
                region=a.sets['plate.vertex_00'],
                u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
            
            self.model.DisplacementBC(name='Pinned-X0Y1', createStepName='Initial', 
                region=a.sets['plate.vertex_01'],
                u1=0.0, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
            
            self.model.DisplacementBC(name='Pinned-X1Y0', createStepName='Initial', 
                region=a.sets['plate.vertex_10'],
                u1=UNSET, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

        #* Apply displacement BCs on the reference points
        self.model.DisplacementBC(name=label_rp, createStepName='Loading', 
            region=a.sets[label_rp],
            u1=self.strain_scale, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)


def extract_field(name_job, fname_save='specimen-field-S4R.dat'):
    
    NAME_INSTANCE = 'PLATE'
    NAME_SET = 'PARTITION_SQUARE'
    
    odb = OdbOperation(name_job)
    element_labels, indices_fieldOutput = odb.get_element_labels_and_indices(
        name_instance=NAME_INSTANCE, name_set=NAME_SET)
    
    coordinates = odb.probe_element_center_coordinate(
        name_instance=NAME_INSTANCE, element_label=element_labels)

    with open(fname_save, 'w') as f:
        
        f.write('Variables= X Y Z index S11 S22 S12 TSHR13 TSHR23 E11 E22 E12\n')
        f.write('zone T=" %s %s "\n'%(NAME_INSTANCE, NAME_SET))

        for i_elem, element_label in enumerate(element_labels):
            
            values_S11 = odb.probe_shell_element_thickness_values(variable='S', component='S11',
                                    name_instance=NAME_INSTANCE, element_label=element_label)
            values_S22 = odb.probe_shell_element_thickness_values(variable='S', component='S22',
                                    name_instance=NAME_INSTANCE, element_label=element_label)
            values_S12 = odb.probe_shell_element_thickness_values(variable='S', component='S12',
                                    name_instance=NAME_INSTANCE, element_label=element_label)
            values_TSHR13 = odb.probe_shell_element_thickness_values(variable='TSHR13', component=None,
                                    name_instance=NAME_INSTANCE, element_label=element_label)
            values_TSHR23 = odb.probe_shell_element_thickness_values(variable='TSHR23', component=None,
                                    name_instance=NAME_INSTANCE, element_label=element_label)
            values_E11 = odb.probe_shell_element_thickness_values(variable='E', component='E11',
                                    name_instance=NAME_INSTANCE, element_label=element_label)
            values_E22 = odb.probe_shell_element_thickness_values(variable='E', component='E22',
                                    name_instance=NAME_INSTANCE, element_label=element_label)
            values_E12 = odb.probe_shell_element_thickness_values(variable='E', component='E12',
                                    name_instance=NAME_INSTANCE, element_label=element_label)

            n_thickness = values_S11.shape[0]
            
            for i_thick in range(1, n_thickness-1):
                f.write(' %14.6E'%(coordinates[i_elem][0]))
                f.write(' %14.6E'%(coordinates[i_elem][1]))
                f.write(' %14.6E'%(coordinates[i_elem][2]+values_S11[i_thick,0]))
                f.write(' %d'%(indices_fieldOutput[i_elem]))
                f.write(' %14.6E'%(values_S11[i_thick,1]))
                f.write(' %14.6E'%(values_S22[i_thick,1]))
                f.write(' %14.6E'%(values_S12[i_thick,1]))
                f.write(' %14.6E'%(values_TSHR13[i_thick-1,1]))
                f.write(' %14.6E'%(values_TSHR23[i_thick-1,1]))
                f.write(' %14.6E'%(values_E11[i_thick,1]))
                f.write(' %14.6E'%(values_E22[i_thick,1]))
                f.write(' %14.6E'%(values_E12[i_thick,1]))
                f.write('\n')

def extract_mid_plane_strain(name_job, fname_save='specimen-mid-plane-strain-S4R.dat'):
    
    NAME_INSTANCE = 'PLATE'
    NAME_SET = 'PARTITION_SQUARE'
    
    odb = OdbOperation(name_job)
    element_labels, indices_fieldOutput = odb.get_element_labels_and_indices(
        name_instance=NAME_INSTANCE, name_set=NAME_SET)
    
    coordinates = odb.probe_element_center_coordinate(
        name_instance=NAME_INSTANCE, element_label=element_labels)

    _element_labels, value_SE = odb.probe_element_set_values(
        step='Loading', frame=-1, variable='SE', component=None,
        name_instance=NAME_INSTANCE, name_set=NAME_SET)
    
    _, value_SF = odb.probe_element_set_values(
        step='Loading', frame=-1, variable='SF', component=None,
        name_instance=NAME_INSTANCE, name_set=NAME_SET)
    
    n_elements = len(element_labels)
    for i in range(n_elements):
        if _element_labels[i] != element_labels[i]:
            raise ValueError('Element labels do not match: %d != %d'%(_element_labels[i], element_labels[i]))

    with open(fname_save, 'w') as f:
        
        f.write('Variables= X Y Z index N11 N22 N12 epsilon11 epsilon22 epsilon12\n')
        f.write('zone T=" %s %s "\n'%(NAME_INSTANCE, NAME_SET))

        for i_elem in range(n_elements):
            for j in range(3):
                f.write(' %14.6E'%(coordinates[i_elem][j]))
            f.write(' %d'%(indices_fieldOutput[i_elem]))
            for j in range(3):
                f.write(' %14.6E'%(value_SF[i_elem, j]))
            for j in range(3):
                f.write(' %14.6E'%(value_SE[i_elem, j]))
            f.write('\n')


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

    name_job = 'Job_OHP_%d_%d'%(index_run, index_strain_vector)

    model = TestModel_InPlaneLoad(name_job, pGeo, pMesh, pRun,
                        strain_component=index_strain_vector,
                        strain_scale=parameters['strain_scale'])
    model.build()
    model.set_view()
    model.save_cae('OHP_%d_%d.cae'%(index_run, index_strain_vector))
    
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
                    f.write('%s_RF  %20.6E \n'%(label_rp, rf_RPs[i_rp]/model.volume_box))
                    
                for i_rp, label_rp in enumerate(model.label_rp):
                    f.write('%s_U   %20.6E \n'%(label_rp, u_RPs[i_rp]))

                f.write('Strain_%d  %20.6E \n'%(index_strain_vector, model.strain_scale))
                f.write('Volume_box  %20.6E \n'%(model.volume_box))
                f.write('Volume_hole %20.6E \n'%(model.volume_hole))
                f.write('Volume      %20.6E \n'%(model.volume))
            
            extract_field(name_job=name_job, fname_save=name_job+'-field.dat')
            extract_mid_plane_strain(name_job=name_job, fname_save=name_job+'-mid-plane.dat')
            