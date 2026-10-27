'''
Classes for applying linear boundary conditions via constraint equations in Abaqus.
'''
import copy

from AbaqusTools import IS_ABAQUS
from AbaqusTools.model import NodeOperation

if IS_ABAQUS:

    from abaqus import *
    from abaqusConstants import *
    from symbolicConstants import *
    from caeModules import *
    import mesh
    import odbAccess


class LBC_3DOrthotropic(NodeOperation):
    '''
    Setup boundary conditions via constraint equations to homogenize
    three-dimensional representative volume elements (RVE) as a 3D orthotropic material.
    
    The special requirement is that the boundary faces are subjected to linear deformations.

    The effective 6x6 stiffness matrix `C` (Voigt notation) is as follows,
    where 1, 2, 3 correspond to x, y, z directions, respectively.

    - {sigma} = C{epsilon}
    - {sigma} = [sigma11, sigma22, sigma33, sigma23, sigma13, sigma12]^T
    - {epsilon} = [epsilon11, epsilon22, epsilon33, gamma23, gamma13, gamma12]^T
    - S = C^{-1} (compliance matrix)
    - gamma_ij = 2 * epsilon_ij (engineering shear strain, i != j)

    Set 6 reference points (RP), use their displacement component `u1` to represent {epsilon}.
    
    - epsilon11 = u1(RP_11)
    - epsilon22 = u1(RP_22)
    - epsilon33 = u1(RP_33)
    - gamma23 = u1(RP_23)
    - gamma13 = u1(RP_13)
    - gamma12 = u1(RP_12)
    
    The engineering constants are:
    
    - E11 = 1 / S11
    - E22 = 1 / S22
    - E33 = 1 / S33
    - G23 = 1 / S44
    - G13 = 1 / S55
    - G12 = 1 / S66
    - niu12 = - S12 / S11
    - niu13 = - S13 / S11
    - niu23 = - S23 / S22

    '''

    @staticmethod
    def create_node_sets(myMdl, name_instance,
                name_master_face_set, name_slave_face_set, coords_sorting,
                name_forbidden_sets=[], label_forbidden_nodes=[],
                name_mfn=None, name_sfn=None):
        return NodeOperation.create_paired_node_sets(myMdl, name_instance, name_master_face_set, name_slave_face_set, coords_sorting, name_forbidden_sets, label_forbidden_nodes, name_mfn, name_sfn)

    @staticmethod
    def create_constraints_strain_vector(myMdl, name_eqn='Lin3D',
                name_mfn_x_set='MFn-X', name_sfn_x_set='SFn-X',
                name_mfn_y_set='MFn-Y', name_sfn_y_set='SFn-Y',
                name_mfn_z_set='MFn-Z', name_sfn_z_set='SFn-Z',
                length_x=1.0, length_y=1.0, length_z=1.0,
                name_rp11='RP_11', name_rp22='RP_22', name_rp33='RP_33',
                name_rp23='RP_23', name_rp13='RP_13', name_rp12='RP_12'):
        '''
        Create constraint equations for specified strain vector in the form of PBC:
        [epsilon11, epsilon22, epsilon33, gamma23, gamma13, gamma12] = 
        [u1(RP_11), u1(RP_22), u1(RP_33), u1(RP_23), u1(RP_13), u1(RP_12)].

        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        name_eqn: str
            name of the equation constraint
            
        name_mfn_x_set, name_mfn_y_set, name_mfn_z_set: str
            name of the set of master face nodes for x, y, z directions
            
        name_sfn_x_set, name_sfn_y_set, name_sfn_z_set: str
            name of the set of slave face nodes for x, y, z directions
            
        length_x, length_y, length_z: float
            the length of the model in x, y, z directions
            
        name_rp11, name_rp22, name_rp33, name_rp23, name_rp13, name_rp12: str
            name of the reference points for the strain vector components

        Notes
        ----------------
        Assuming small deformation, one corner of the RVE is at the origin,
        and the faces goes through the origin are called the slave faces,
        the boundary conditions are:
        
        - u1(Mx) = u1(RP_11) * length_x + 0.5*( u1(RP_12) * y(Mx)    + u1(RP_13) * z(Mx) )
        - u2(Mx) = u1(RP_22) * y(Mx)    + 0.5*( u1(RP_12) * length_x + u1(RP_23) * z(Mx) )
        - u3(Mx) = u1(RP_33) * z(Mx)    + 0.5*( u1(RP_13) * length_x + u1(RP_23) * y(Mx) )
        
        - u1(Sx) = 0                    + 0.5*( u1(RP_12) * y(Sx)    + u1(RP_13) * z(Sx) )
        - u2(Sx) = u1(RP_22) * y(Sx)    + 0.5*( 0                    + u1(RP_23) * z(Sx) )
        - u3(Sx) = u1(RP_33) * z(Sx)    + 0.5*( 0                    + u1(RP_23) * y(Sx) )

        - u1(My) = u1(RP_11) * x(My)    + 0.5*( u1(RP_12) * length_y + u1(RP_13) * z(My) )
        - u2(My) = u1(RP_22) * length_y + 0.5*( u1(RP_12) * x(My)    + u1(RP_23) * z(My) )
        - u3(My) = u1(RP_33) * z(My)    + 0.5*( u1(RP_23) * length_y + u1(RP_13) * x(My) )

        - u1(Sy) = u1(RP_11) * x(Sy)    + 0.5*( 0                    + u1(RP_13) * z(Sy) )
        - u2(Sy) = 0                    + 0.5*( u1(RP_12) * x(Sy)    + u1(RP_23) * z(Sy) )
        - u3(Sy) = u1(RP_33) * z(Sy)    + 0.5*( u1(RP_13) * x(Sy)    + 0                 )

        - u1(Mz) = u1(RP_11) * x(Mz)    + 0.5*( u1(RP_13) * length_z + u1(RP_12) * y(Mz) )
        - u2(Mz) = u1(RP_22) * y(Mz)    + 0.5*( u1(RP_23) * length_z + u1(RP_12) * x(Mz) )
        - u3(Mz) = u1(RP_33) * length_z + 0.5*( u1(RP_13) * x(Mz)    + u1(RP_23) * y(Mz) )
        
        - u1(Sz) = u1(RP_11) * x(Sz)    + 0.5*( u1(RP_12) * y(Sz)    + 0                 )
        - u2(Sz) = u1(RP_22) * y(Sz)    + 0.5*( u1(RP_12) * x(Sz)    + 0                 )
        - u3(Sz) = 0                    + 0.5*( u1(RP_13) * x(Sz)    + u1(RP_23) * y(Sz) )
        '''
        aa = myMdl.rootAssembly
        
        #* X-direction: master face
        
        face_nodes = aa.sets[name_mfn_x_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_M = '%s-%d'%(name_mfn_x_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)

            #* Coordinates
            y = face_nodes[i_node].coordinates[1]
            z = face_nodes[i_node].coordinates[2]
            
            #* Constraint equations
            '''
            - u1(Mx) = u1(RP_11) * length_x + 0.5*( u1(RP_12) * y(Mx)    + u1(RP_13) * z(Mx) )
            - u2(Mx) = u2(RP_22) * y(Mx)    + 0.5*( u1(RP_12) * length_x + u1(RP_23) * z(Mx) )
            - u3(Mx) = u3(RP_33) * z(Mx)    + 0.5*( u1(RP_13) * length_x + u1(RP_23) * y(Mx) )
            '''
            myMdl.Equation(name='%s_MX-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 1), (-length_x, name_rp11, 1), (-0.5*y, name_rp12, 1), (-0.5*z, name_rp13, 1)))
            
            myMdl.Equation(name='%s_MX-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 2), (-y, name_rp22, 1), (-0.5*length_x, name_rp12, 1), (-0.5*z, name_rp23, 1)))
            
            myMdl.Equation(name='%s_MX-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 3), (-z, name_rp33, 1), (-0.5*length_x, name_rp13, 1), (-0.5*y, name_rp23, 1)))
    
        #* X-direction: slave face
        
        face_nodes = aa.sets[name_sfn_x_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_S = '%s-%d'%(name_sfn_x_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_S)

            #* Coordinates
            y = face_nodes[i_node].coordinates[1]
            z = face_nodes[i_node].coordinates[2]
            
            #* Constraint equations
            '''
            - u1(Sx) = 0                    + 0.5*( u1(RP_12) * y(Sx)    + u1(RP_13) * z(Sx) )
            - u2(Sx) = u2(RP_22) * y(Sx)    + 0.5*( 0                    + u1(RP_23) * z(Sx) )
            - u3(Sx) = u3(RP_33) * z(Sx)    + 0.5*( 0                    + u1(RP_23) * y(Sx) )
            '''
            myMdl.Equation(name='%s_SX-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 1), (-0.5*y, name_rp12, 1), (-0.5*z, name_rp13, 1)))
            
            myMdl.Equation(name='%s_SX-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 2), (-y, name_rp22, 1), (-0.5*z, name_rp23, 1)))
            
            myMdl.Equation(name='%s_SX-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 3), (-z, name_rp33, 1), (-0.5*y, name_rp23, 1)))

        #* Y-direction: master face
        
        face_nodes = aa.sets[name_mfn_y_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_M = '%s-%d'%(name_mfn_y_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)

            #* Coordinates
            x = face_nodes[i_node].coordinates[0]
            z = face_nodes[i_node].coordinates[2]

            #* Constraint equations
            '''
            - u1(My) = u1(RP_11) * x(My)    + 0.5*( u1(RP_12) * length_y + u1(RP_13) * z(My) )
            - u2(My) = u2(RP_22) * length_y + 0.5*( u1(RP_12) * x(My)    + u1(RP_23) * z(My) )
            - u3(My) = u3(RP_33) * z(My)    + 0.5*( u1(RP_23) * length_y + u1(RP_13) * x(My) )
            '''
            myMdl.Equation(name='%s_MY-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 1), (-x, name_rp11, 1), (-0.5*length_y, name_rp12, 1), (-0.5*z, name_rp13, 1)))
            
            myMdl.Equation(name='%s_MY-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 2), (-length_y, name_rp22, 1), (-0.5*x, name_rp12, 1), (-0.5*z, name_rp23, 1)))
            
            myMdl.Equation(name='%s_MY-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 3), (-z, name_rp33, 1), (-0.5*length_y, name_rp23, 1), (-0.5*x, name_rp13, 1)))
    
        #* Y-direction: slave face
        
        face_nodes = aa.sets[name_sfn_y_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_S = '%s-%d'%(name_sfn_y_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_S)

            #* Coordinates
            x = face_nodes[i_node].coordinates[0]
            z = face_nodes[i_node].coordinates[2]
            
            #* Constraint equations
            '''
            - u1(Sy) = u1(RP_11) * x(Sy)    + 0.5*( 0                    + u1(RP_13) * z(Sy) )
            - u2(Sy) = 0                    + 0.5*( u1(RP_12) * x(My)    + u1(RP_23) * z(Sy) )
            - u3(Sy) = u3(RP_33) * z(Sy)    + 0.5*( u1(RP_13) * x(Sy)    + 0                 )
            '''
            myMdl.Equation(name='%s_SY-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 1), (-x, name_rp11, 1), (-0.5*z, name_rp13, 1)))
            
            myMdl.Equation(name='%s_SY-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 2), (-0.5*x, name_rp12, 1), (-0.5*z, name_rp23, 1)))
            
            myMdl.Equation(name='%s_SY-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 3), (-z, name_rp33, 1), (-0.5*x, name_rp13, 1)))

        #* Z-direction: master face
        
        face_nodes = aa.sets[name_mfn_z_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_M = '%s-%d'%(name_mfn_z_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)
            
            #* Coordinates
            x = face_nodes[i_node].coordinates[0]
            y = face_nodes[i_node].coordinates[1]
            
            #* Constraint equations
            '''
            - u1(Mz) = u1(RP_11) * x(Mz)    + 0.5*( u1(RP_13) * length_z + u1(RP_12) * y(Mz) )
            - u2(Mz) = u2(RP_22) * y(Mz)    + 0.5*( u1(RP_23) * length_z + u1(RP_12) * x(Mz) )
            - u3(Mz) = u3(RP_33) * length_z + 0.5*( u1(RP_13) * x(Mz)    + u1(RP_23) * y(Mz) )
            '''
            myMdl.Equation(name='%s_MZ-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 1), (-x, name_rp11, 1), (-0.5*length_z, name_rp13, 1), (-0.5*y, name_rp12, 1)))
            
            myMdl.Equation(name='%s_MZ-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 2), (-y, name_rp22, 1), (-0.5*length_z, name_rp23, 1), (-0.5*x, name_rp12, 1)))
            
            myMdl.Equation(name='%s_MZ-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 3), (-length_z, name_rp33, 1), (-0.5*x, name_rp13, 1), (-0.5*y, name_rp23, 1)))

        #* Z-direction: slave face
        
        face_nodes = aa.sets[name_sfn_z_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_S = '%s-%d'%(name_sfn_z_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_S)
            
            #* Coordinates
            x = face_nodes[i_node].coordinates[0]
            y = face_nodes[i_node].coordinates[1]
            
            #* Constraint equations
            '''
            - u1(Sz) = u1(RP_11) * x(Sz)    + 0.5*( u1(RP_12) * y(Sz)    + 0                 )
            - u2(Sz) = u2(RP_22) * y(Sz)    + 0.5*( u1(RP_12) * x(Sz)    + 0                 )
            - u3(Sz) = 0                    + 0.5*( u1(RP_13) * x(Sz)    + u1(RP_23) * y(Sz) )
            '''
            myMdl.Equation(name='%s_SZ-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 1), (-x, name_rp11, 1), (-0.5*y, name_rp12, 1)))
            
            myMdl.Equation(name='%s_SZ-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 2), (-y, name_rp22, 1), (-0.5*x, name_rp12, 1)))
            
            myMdl.Equation(name='%s_SZ-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 3), (-0.5*x, name_rp13, 1), (-0.5*y, name_rp23, 1)))
        
    @staticmethod
    def calculate_engineering_constants(stiffness_matrix):
        '''
        Calculate the engineering constants from the stiffness matrix.
        
        Parameters
        ---------------
        stiffness_matrix: ndarray [6, 6]
            the stiffness matrix
            
        Returns
        ---------------
        engineering_constants: dict
            the engineering constants
            - E11: float
            - E22: float
            - E33: float
            - G23: float
            - G13: float
            - G12: float
            - niu12: float
            - niu13: float
            - niu23: float
            - C_avg: list [6, 6]
                the average stiffness matrix (made to be symmetric)
            - compliance_matrix: list [6, 6]
                the compliance matrix
            - stiffness_matrix: list [6, 6]
                the original stiffness matrix
        '''
        import numpy as np
        
        stiffness_matrix = np.array(stiffness_matrix)
        C_avg = 0.5*(stiffness_matrix + stiffness_matrix.T)
        S_avg = np.linalg.inv(C_avg)
        
        result = {
            'E11': 1.0 / S_avg[0, 0],
            'E22': 1.0 / S_avg[1, 1],
            'E33': 1.0 / S_avg[2, 2],
            'G23': 1.0 / S_avg[3, 3],
            'G13': 1.0 / S_avg[4, 4],
            'G12': 1.0 / S_avg[5, 5],
            'niu12': - S_avg[0, 1] / S_avg[0, 0],
            'niu13': - S_avg[0, 2] / S_avg[0, 0],
            'niu23': - S_avg[1, 2] / S_avg[1, 1],
            'stiffness_matrix': stiffness_matrix.tolist(),
            'compliance_matrix': S_avg.tolist(),
            'C_avg': C_avg.tolist()
        }
        return result
    
    @staticmethod
    def calculate_rp_displacements(strain_vector, model_length_x, model_length_y, model_length_z):
        '''
        Calculate the displacement `u` of reference points for boundary condition definition.

        Parameters
        ----------------------
        strain_vector: ndarray [6]
            the strain vector for loading
            
        model_length_x, model_length_y, model_length_z: float
            the length of model in x, y, z directions
        
        Returns
        ----------------------
        displacements: dict
            the `u1` displacement of reference points
            - RP_11: float
            - RP_22: float
            - RP_33: float
            - RP_23: float
            - RP_13: float
            - RP_12: float
        '''
        volume_box = model_length_x*model_length_y*model_length_z
            
        displacements = {
            'RP_11': strain_vector[0]*volume_box,
            'RP_22': strain_vector[1]*volume_box,
            'RP_33': strain_vector[2]*volume_box,
            'RP_23': strain_vector[3]*volume_box,
            'RP_13': strain_vector[4]*volume_box,
            'RP_12': strain_vector[5]*volume_box,
        }

        return displacements


class LBC_3DOrthotropic_2(LBC_3DOrthotropic):
    '''
    
    '''

    @staticmethod
    def exclude_forbidden_nodes_pbc(myMdl, name_instance,
                        master_face_nodes, slave_face_nodes,
                        name_forbidden_sets=[], label_forbidden_nodes=[]):
        '''
        Exclude forbidden nodes from the master face,
        and exclude the corresponding node from the slave face.
        
        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        name_instance: str
            name of the instance that the set belongs to.
            
        master_face_nodes: MeshNodeArray object
            nodes on the master face
            
        slave_face_nodes: MeshNodeArray object
            nodes on the slave face
            
        name_forbidden_sets: list of strings
            name of sets that consists of forbidden nodes.
            It can be a node set, or a geometry set of vertices, edges, or faces.
            
        label_forbidden_nodes: list of integers
            labels of forbidden nodes.
            An overlap between `name_forbidden_sets` and `label_forbidden_nodes` is allowed.
            
        Return
        ---------------
        label_master_nodes: list of integers
            labels of the master nodes
            
        label_slave_nodes: list of integers
            labels of the slave nodes
            
        label_forbidden: list of integers
            labels of forbidden nodes
            
        Note
        ---------------
        The forbidden nodes only contain nodes in the master faces.
        
        Examples
        ---------------
        >>> label_master_nodes, label_slave_nodes, label_forbidden = \\
        >>>         PeriodicBC.exclude_forbidden_nodes_pbc(myMdl, name_instance, master_face_nodes, 
        >>>         slave_face_nodes, name_forbidden_sets=[], label_forbidden_nodes=[])
        '''
        label_forbidden = copy.deepcopy(label_forbidden_nodes)
        label_master_nodes = []
        label_slave_nodes = []

        for name_set in name_forbidden_sets:
            
            set_ = myMdl.rootAssembly.sets['%s.%s'%(name_instance, name_set)]
            label_forbidden += [node.label for node in set_.nodes]

        for mfn, sfn in zip(master_face_nodes, slave_face_nodes):

            # Check master node is not among the forbidden nodes list
            if mfn.label in label_forbidden:
                continue
            
            label_master_nodes.append(mfn.label)
            label_slave_nodes.append(sfn.label)
            
            # Update forbidden nodes list
            label_forbidden.append(mfn.label)
                
        return label_master_nodes, label_slave_nodes, label_forbidden

    @staticmethod
    def create_node_sets(myMdl, name_instance,
                        name_master_face_set, name_slave_face_set, coords_sorting,
                        name_forbidden_sets=[], label_forbidden_nodes=[],
                        name_mfn=None, name_sfn=None, bc_type='normal'):
        '''
        Create node sets on the master/slave faces for the linear boundary condition.
        
        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        name_instance: str
            name of the instance that the set belongs to.
            
        name_master_face_set: str
            name of the master face set
            
        name_slave_face_set: str
            name of the slave face set
            
        coords_sorting: tuple of int
            a tuple contains coordinate indices for sorting nodes, e.g., (0,1), (1,2), (2,0).
            
        name_forbidden_sets: list of strings
            name of sets that consists of forbidden nodes.
            It can be a node set, or a geometry set of vertices, edges, or faces.
            
        label_forbidden_nodes: list of integers
            labels of forbidden nodes.
            An overlap between `name_forbidden_sets` and `label_forbidden_nodes` is allowed.
            
        name_mfn, name_sfn: str
            name of the node sets on the master/slave faces.
            If None, use default name.
            
        bc_type: str
            type of boundary condition, 'normal' or 'shear'.
            
        Return
        ----------------
        name_mfn, name_sfn: str
            name of the node sets on the master/slave faces.
        '''
        #* Sort nodes by coordinates to find node pairs automatically
        master_face_nodes = LBC_3DOrthotropic.get_nodes_from_face(myMdl, name_master_face_set, coords_sorting, name_instance)
        slave_face_nodes  = LBC_3DOrthotropic.get_nodes_from_face(myMdl, name_slave_face_set,  coords_sorting, name_instance)
        
        num_node_master = len(master_face_nodes)
        num_node_slave  = len(slave_face_nodes)
                
        #* Exclude forbidden nodes
        if bc_type in ['normal', 'pbc_z'] and 'z' in name_master_face_set:
        
            labels_master, labels_slave, label_forbidden_nodes = \
                LBC_3DOrthotropic_2.exclude_forbidden_nodes_pbc(
                    myMdl, name_instance, master_face_nodes, slave_face_nodes,
                    name_forbidden_sets, label_forbidden_nodes)
                
        else:
            
            labels_master, label_forbidden_nodes = \
                LBC_3DOrthotropic.exclude_forbidden_nodes(myMdl, nodes=master_face_nodes, 
                                    name_forbidden_sets=name_forbidden_sets, 
                                    label_forbidden_nodes=label_forbidden_nodes, 
                                    name_part=None, name_instance=name_instance)
                
            labels_slave, label_forbidden_nodes = \
                LBC_3DOrthotropic.exclude_forbidden_nodes(myMdl, nodes=slave_face_nodes, 
                                    name_forbidden_sets=name_forbidden_sets, 
                                    label_forbidden_nodes=label_forbidden_nodes, 
                                    name_part=None, name_instance=name_instance)
        
        #* Node set names
        if name_mfn is None:
            name_mfn = 'MFNode-%s-%s'%(name_instance, name_master_face_set)
        if name_sfn is None:
            name_sfn = 'SFNode-%s-%s'%(name_instance, name_slave_face_set)
                
        #* Create node sets
        myMdl.rootAssembly.SetFromNodeLabels(name=name_mfn, nodeLabels=((name_instance, labels_master),), unsorted=True)
        myMdl.rootAssembly.SetFromNodeLabels(name=name_sfn, nodeLabels=((name_instance, labels_slave ),), unsorted=True)
        
        #* Print information
        print('>>> --------------------')
        print('[Linear boundary conditions] Instance: %s'%(name_instance))
        print('    Node sets: master face [%s] <=> slave face [%s];'%(name_master_face_set, name_slave_face_set))
        print('    Master face: nNodes= %d after excluding forbidden nodes from %d nodes'%(len(labels_master), num_node_master))
        print('    Slave face:  nNodes= %d after excluding forbidden nodes from %d nodes'%(len(labels_slave), num_node_slave))
        print('>>>')
        
        return name_mfn, name_sfn, label_forbidden_nodes

    @staticmethod
    def create_constraints_strain_vector(myMdl, name_eqn='Lin3D_2',
                name_mfn_x_set='MFn-X', name_sfn_x_set='SFn-X',
                name_mfn_y_set='MFn-Y', name_sfn_y_set='SFn-Y',
                name_mfn_z_set='MFn-Z', name_sfn_z_set='SFn-Z',
                length_x=1.0, length_y=1.0, length_z=1.0,
                name_rp11='RP_11', name_rp22='RP_22', name_rp33='RP_33',
                name_rp23='RP_23', name_rp13='RP_13', name_rp12='RP_12',
                bc_type='normal'):
        '''
        Create constraint equations for specified strain vector in the form of PBC:
        [epsilon11, epsilon22, epsilon33, gamma23, gamma13, gamma12] = 
        [u1(RP_11), u1(RP_22), u1(RP_33), u1(RP_23), u1(RP_13), u1(RP_12)].

        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        name_eqn: str
            name of the equation constraint
            
        name_mfn_x_set, name_mfn_y_set, name_mfn_z_set: str
            name of the set of master face nodes for x, y, z directions
            
        name_sfn_x_set, name_sfn_y_set, name_sfn_z_set: str
            name of the set of slave face nodes for x, y, z directions
            
        length_x, length_y, length_z: float
            the length of the model in x, y, z directions
            
        name_rp11, name_rp22, name_rp33, name_rp23, name_rp13, name_rp12: str
            name of the reference points for the strain vector components
            
        bc_type: str
            type of boundary condition, 'normal' or 'shear'.

        Notes
        ----------------
        Assuming small deformation, one corner of the RVE is at the origin,
        and the faces goes through the origin are called the slave faces,
        the boundary conditions are:
        
        - u1(Mx) = u1(RP_11) * length_x + 0.5*( u1(RP_12) * y(Mx)    + u1(RP_13) * z(Mx) )
        - u2(Mx) = u2(RP_22) * y(Mx)    + 0.5*( u1(RP_12) * length_x + u1(RP_23) * z(Mx) )
        - u3(Mx) = u3(RP_33) * z(Mx)    + 0.5*( u1(RP_13) * length_x + u1(RP_23) * y(Mx) )
        
        - u1(Sx) = 0                    + 0.5*( u1(RP_12) * y(Sx)    + u1(RP_13) * z(Sx) )
        - u2(Sx) = u2(RP_22) * y(Sx)    + 0.5*( 0                    + u1(RP_23) * z(Sx) )
        - u3(Sx) = u3(RP_33) * z(Sx)    + 0.5*( 0                    + u1(RP_23) * y(Sx) )

        - u1(My) = u1(RP_11) * x(My)    + 0.5*( u1(RP_12) * length_y + u1(RP_13) * z(My) )
        - u2(My) = u2(RP_22) * length_y + 0.5*( u1(RP_12) * x(My)    + u1(RP_23) * z(My) )
        - u3(My) = u3(RP_33) * z(My)    + 0.5*( u1(RP_23) * length_y + u1(RP_13) * x(My) )

        - u1(Sy) = u1(RP_11) * x(Sy)    + 0.5*( 0                    + u1(RP_13) * z(Sy) )
        - u2(Sy) = 0                    + 0.5*( u1(RP_12) * x(My)    + u1(RP_23) * z(Sy) )
        - u3(Sy) = u3(RP_33) * z(Sy)    + 0.5*( u1(RP_13) * x(Sy)    + 0                 )

        Z-direction faces subject to periodic boundary condition,
        when under normal strain state (epsilon_ii):
        
        - u1(Mz) - u1(Sz) = 0
        - u2(Mz) - u2(Sz) = 0
        - u3(Mz) - u3(Sz) - length_z * u1(RP_33) = 0

        Z-direction faces subject to linear boundary condition,
        when under shear strain state (gamma_ij):
        
        - u1(Mz) = 0.5*( u1(RP_13) * length_z + u1(RP_12) * y(Mz) )
        - u2(Mz) = 0.5*( u1(RP_23) * length_z + u1(RP_12) * x(Mz) )
        - u3(Mz) = 0.5*( u1(RP_13) * x(Mz)    + u1(RP_23) * y(Mz) )
        
        - u1(Sz) = 0.5*( u1(RP_12) * y(Sz)    + 0                 )
        - u2(Sz) = 0.5*( u1(RP_12) * x(Sz)    + 0                 )
        - u3(Sz) = 0.5*( u1(RP_13) * x(Sz)    + u1(RP_23) * y(Sz) )
        '''
        aa = myMdl.rootAssembly
        
        #* X-direction: master face
        
        face_nodes = aa.sets[name_mfn_x_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_M = '%s-%d'%(name_mfn_x_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)

            #* Coordinates
            y = face_nodes[i_node].coordinates[1]
            z = face_nodes[i_node].coordinates[2]
            
            #* Constraint equations
            '''
            - u1(Mx) = u1(RP_11) * length_x + 0.5*( u1(RP_12) * y(Mx)    + u1(RP_13) * z(Mx) )
            - u2(Mx) = u2(RP_22) * y(Mx)    + 0.5*( u1(RP_12) * length_x + u1(RP_23) * z(Mx) )
            - u3(Mx) = u3(RP_33) * z(Mx)    + 0.5*( u1(RP_13) * length_x + u1(RP_23) * y(Mx) )
            '''
            myMdl.Equation(name='%s_MX-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 1), (-length_x, name_rp11, 1), (-0.5*y, name_rp12, 1), (-0.5*z, name_rp13, 1)))
            
            myMdl.Equation(name='%s_MX-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 2), (-y, name_rp22, 1), (-0.5*length_x, name_rp12, 1), (-0.5*z, name_rp23, 1)))
            
            myMdl.Equation(name='%s_MX-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 3), (-z, name_rp33, 1), (-0.5*length_x, name_rp13, 1), (-0.5*y, name_rp23, 1)))
    
        #* X-direction: slave face
        
        face_nodes = aa.sets[name_sfn_x_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_S = '%s-%d'%(name_sfn_x_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_S)

            #* Coordinates
            y = face_nodes[i_node].coordinates[1]
            z = face_nodes[i_node].coordinates[2]
            
            #* Constraint equations
            '''
            - u1(Sx) = 0                    + 0.5*( u1(RP_12) * y(Sx)    + u1(RP_13) * z(Sx) )
            - u2(Sx) = u2(RP_22) * y(Sx)    + 0.5*( 0                    + u1(RP_23) * z(Sx) )
            - u3(Sx) = u3(RP_33) * z(Sx)    + 0.5*( 0                    + u1(RP_23) * y(Sx) )
            '''
            myMdl.Equation(name='%s_SX-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 1), (-0.5*y, name_rp12, 1), (-0.5*z, name_rp13, 1)))
            
            myMdl.Equation(name='%s_SX-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 2), (-y, name_rp22, 1), (-0.5*z, name_rp23, 1)))
            
            myMdl.Equation(name='%s_SX-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 3), (-z, name_rp33, 1), (-0.5*y, name_rp23, 1)))

        #* Y-direction: master face
        
        face_nodes = aa.sets[name_mfn_y_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_M = '%s-%d'%(name_mfn_y_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)

            #* Coordinates
            x = face_nodes[i_node].coordinates[0]
            z = face_nodes[i_node].coordinates[2]

            #* Constraint equations
            '''
            - u1(My) = u1(RP_11) * x(My)    + 0.5*( u1(RP_12) * length_y + u1(RP_13) * z(My) )
            - u2(My) = u2(RP_22) * length_y + 0.5*( u1(RP_12) * x(My)    + u1(RP_23) * z(My) )
            - u3(My) = u3(RP_33) * z(My)    + 0.5*( u1(RP_23) * length_y + u1(RP_13) * x(My) )
            '''
            myMdl.Equation(name='%s_MY-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 1), (-x, name_rp11, 1), (-0.5*length_y, name_rp12, 1), (-0.5*z, name_rp13, 1)))
            
            myMdl.Equation(name='%s_MY-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 2), (-length_y, name_rp22, 1), (-0.5*x, name_rp12, 1), (-0.5*z, name_rp23, 1)))
            
            myMdl.Equation(name='%s_MY-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 3), (-z, name_rp33, 1), (-0.5*length_y, name_rp23, 1), (-0.5*x, name_rp13, 1)))
    
        #* Y-direction: slave face
        
        face_nodes = aa.sets[name_sfn_y_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_S = '%s-%d'%(name_sfn_y_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_S)

            #* Coordinates
            x = face_nodes[i_node].coordinates[0]
            z = face_nodes[i_node].coordinates[2]
            
            #* Constraint equations
            '''
            - u1(Sy) = u1(RP_11) * x(Sy)    + 0.5*( 0                    + u1(RP_13) * z(Sy) )
            - u2(Sy) = 0                    + 0.5*( u1(RP_12) * x(My)    + u1(RP_23) * z(Sy) )
            - u3(Sy) = u3(RP_33) * z(Sy)    + 0.5*( u1(RP_13) * x(Sy)    + 0                 )
            '''
            myMdl.Equation(name='%s_SY-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 1), (-x, name_rp11, 1), (-0.5*z, name_rp13, 1)))
            
            myMdl.Equation(name='%s_SY-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 2), (-0.5*x, name_rp12, 1), (-0.5*z, name_rp23, 1)))
            
            myMdl.Equation(name='%s_SY-%d-3'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 3), (-z, name_rp33, 1), (-0.5*x, name_rp13, 1)))

        #* Z-direction

        if bc_type == 'normal':

            master_face_nodes = aa.sets[name_mfn_z_set].nodes
            slave_face_nodes  = aa.sets[name_sfn_z_set].nodes
            
            for i_node in range(len(master_face_nodes)):

                #* Create sets for each pair of nodes
                name_M = '%s-%d'%(name_mfn_z_set, i_node)
                name_S = '%s-%d'%(name_sfn_z_set, i_node)
                aa.Set(nodes=mesh.MeshNodeArray((master_face_nodes[i_node],)), name=name_M)
                aa.Set(nodes=mesh.MeshNodeArray((slave_face_nodes [i_node],)), name=name_S)

                #* Constraint equations
                '''
                - u1(Mz) - u1(Sz) = 0
                - u2(Mz) - u2(Sz) = 0
                - u3(Mz) - u3(Sz) - length_z * u1(RP_33) = 0
                '''
                myMdl.Equation(name='%s_Z-%d-1'%(name_eqn, i_node), 
                        terms=((1.0, name_M, 1), (-1.0, name_S, 1)))
                
                myMdl.Equation(name='%s_Z-%d-2'%(name_eqn, i_node), 
                        terms=((1.0, name_M, 2), (-1.0, name_S, 2)))
                
                myMdl.Equation(name='%s_Z-%d-3'%(name_eqn, i_node), 
                        terms=((1.0, name_M, 3), (-1.0, name_S, 3), (-length_z, name_rp33, 1)))

        elif bc_type == 'shear':

            #* Z-direction: master face
            
            face_nodes = aa.sets[name_mfn_z_set].nodes

            for i_node in range(len(face_nodes)):

                #* Create set for node
                name_M = '%s-%d'%(name_mfn_z_set, i_node)
                aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)
                
                #* Coordinates
                x = face_nodes[i_node].coordinates[0]
                y = face_nodes[i_node].coordinates[1]
                
                #* Constraint equations
                '''
                - u1(Mz) = 0.5*( u1(RP_13) * length_z + u1(RP_12) * y(Mz) )
                - u2(Mz) = 0.5*( u1(RP_23) * length_z + u1(RP_12) * x(Mz) )
                - u3(Mz) = 0.5*( u1(RP_13) * x(Mz)    + u1(RP_23) * y(Mz) )
                '''
                myMdl.Equation(name='%s_MZ-%d-1'%(name_eqn, i_node), 
                        terms=((1.0, name_M, 1), (-0.5*length_z, name_rp13, 1), (-0.5*y, name_rp12, 1)))
                
                myMdl.Equation(name='%s_MZ-%d-2'%(name_eqn, i_node), 
                        terms=((1.0, name_M, 2), (-0.5*length_z, name_rp23, 1), (-0.5*x, name_rp12, 1)))
                
                myMdl.Equation(name='%s_MZ-%d-3'%(name_eqn, i_node), 
                        terms=((1.0, name_M, 3), (-0.5*x, name_rp13, 1), (-0.5*y, name_rp23, 1)))

            #* Z-direction: slave face
            
            face_nodes = aa.sets[name_sfn_z_set].nodes

            for i_node in range(len(face_nodes)):

                #* Create set for node
                name_S = '%s-%d'%(name_sfn_z_set, i_node)
                aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_S)
                
                #* Coordinates
                x = face_nodes[i_node].coordinates[0]
                y = face_nodes[i_node].coordinates[1]
                
                #* Constraint equations
                '''
                - u1(Sz) = 0.5*( u1(RP_12) * y(Sz)    + 0                 )
                - u2(Sz) = 0.5*( u1(RP_12) * x(Sz)    + 0                 )
                - u3(Sz) = 0.5*( u1(RP_13) * x(Sz)    + u1(RP_23) * y(Sz) )
                '''
                myMdl.Equation(name='%s_SZ-%d-1'%(name_eqn, i_node), 
                        terms=((1.0, name_S, 1), (-0.5*y, name_rp12, 1)))
                
                myMdl.Equation(name='%s_SZ-%d-2'%(name_eqn, i_node), 
                        terms=((1.0, name_S, 2), (-0.5*x, name_rp12, 1)))
                
                myMdl.Equation(name='%s_SZ-%d-3'%(name_eqn, i_node), 
                        terms=((1.0, name_S, 3), (-0.5*x, name_rp13, 1), (-0.5*y, name_rp23, 1)))
            
        elif bc_type == 'pbc_z':
            
            master_face_nodes = aa.sets[name_mfn_z_set].nodes
            slave_face_nodes  = aa.sets[name_sfn_z_set].nodes
            
            for i_node in range(len(master_face_nodes)):

                #* Create sets for each pair of nodes
                name_M = '%s-%d'%(name_mfn_z_set, i_node)
                name_S = '%s-%d'%(name_sfn_z_set, i_node)
                aa.Set(nodes=mesh.MeshNodeArray((master_face_nodes[i_node],)), name=name_M)
                aa.Set(nodes=mesh.MeshNodeArray((slave_face_nodes [i_node],)), name=name_S)

                #* Constraint equations
                '''
                - u1(Mz) - u1(Sz) - 0.5 * length_z * u1(RP_13) = 0
                - u2(Mz) - u2(Sz) - 0.5 * length_z * u1(RP_23) = 0
                - u3(Mz) - u3(Sz) - length_z * u1(RP_33) = 0
                '''
                myMdl.Equation(name='%s_Z-%d-1'%(name_eqn, i_node), 
                        terms=((1.0, name_M, 1), (-1.0, name_S, 1), (-0.5*length_z, name_rp13, 1)))
                
                myMdl.Equation(name='%s_Z-%d-2'%(name_eqn, i_node), 
                        terms=((1.0, name_M, 2), (-1.0, name_S, 2), (-0.5*length_z, name_rp23, 1)))
                
                myMdl.Equation(name='%s_Z-%d-3'%(name_eqn, i_node), 
                        terms=((1.0, name_M, 3), (-1.0, name_S, 3), (-length_z, name_rp33, 1)))
            
        else:
            raise ValueError("Invalid boundary condition type: %s"%(bc_type))


class LBC_InPlaneLoad(NodeOperation):
    '''
    Linear boundary condition for in-plane loading of 3D models.
    The out-of-plane faces are free of constraints, and the in-plane faces are subject to linear boundary conditions.
    
    Set 3 reference points (RP), use their displacement component `u1` to represent {epsilon}.
    
    - epsilon11 = u1(RP_11)
    - epsilon22 = u1(RP_22)
    - gamma12 = u1(RP_12)
    
    ''' 
    @staticmethod
    def create_node_sets(myMdl, name_instance,
                name_master_face_set, name_slave_face_set, coords_sorting,
                name_forbidden_sets=[], label_forbidden_nodes=[],
                name_mfn=None, name_sfn=None):
        return NodeOperation.create_paired_node_sets(myMdl, name_instance,
                name_master_face_set, name_slave_face_set, coords_sorting,
                name_forbidden_sets, label_forbidden_nodes, name_mfn, name_sfn)
    
    @staticmethod
    def create_constraints_gamma12(myMdl, name_eqn='InPlane12',
                name_mfn_x_set='MFn-X', name_sfn_x_set='SFn-X',
                name_mfn_y_set='MFn-Y', name_sfn_y_set='SFn-Y',
                length_x=1.0, length_y=1.0,
                name_rp12='RP_12'):
        '''
        Create constraint equations for gamma12 = u1(RP_12).

        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        name_eqn: str
            name of the equation constraint
            
        name_mfn_x_set, name_mfn_y_set: str
            name of the set of master face nodes for x, y directions
            
        name_sfn_x_set, name_sfn_y_set: str
            name of the set of slave face nodes for x, y directions
            
        length_x, length_y: float
            the length of the model in x, y directions
            
        name_rp12: str
            name of the reference point for the strain vector component
            
        Notes
        ----------------
        Assuming small deformation, one corner of the RVE is at the origin,
        and the faces goes through the origin are called the slave faces,
        the boundary conditions are:
        
        - u1(Mx) = 0.5*( u1(RP_12) * y(Mx)    )
        - u2(Mx) = 0.5*( u1(RP_12) * length_x )
        
        - u1(Sx) = 0.5*( u1(RP_12) * y(Sx)    )
        - u2(Sx) = 0

        - u1(My) = 0.5*( u1(RP_12) * length_y )
        - u2(My) = 0.5*( u1(RP_12) * x(My)    )

        - u1(Sy) = 0
        - u2(Sy) = 0.5*( u1(RP_12) * x(Sy)    )
        '''
        aa = myMdl.rootAssembly
        
        #* X-direction: master face
        
        face_nodes = aa.sets[name_mfn_x_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_M = '%s-%d'%(name_mfn_x_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)

            #* Coordinates
            y = face_nodes[i_node].coordinates[1]

            #* Constraint equations
            '''
            - u1(Mx) = 0.5*( u1(RP_12) * y(Mx)    )
            - u2(Mx) = 0.5*( u1(RP_12) * length_x )
            '''
            myMdl.Equation(name='%s_MX-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 1), (-0.5*y, name_rp12, 1)))
            
            myMdl.Equation(name='%s_MX-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 2), (-0.5*length_x, name_rp12, 1)))
    
        #* X-direction: slave face
        
        face_nodes = aa.sets[name_sfn_x_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_S = '%s-%d'%(name_sfn_x_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_S)

            #* Coordinates
            y = face_nodes[i_node].coordinates[1]
            
            #* Constraint equations
            '''
            - u1(Sx) = 0.5*( u1(RP_12) * y(Sx)    )
            - u2(Sx) = 0
            '''
            myMdl.Equation(name='%s_SX-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 1), (-0.5*y, name_rp12, 1)))
            
            myMdl.Equation(name='%s_SX-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 2), (0.0, name_rp12, 1)))

        #* Y-direction: master face
        
        face_nodes = aa.sets[name_mfn_y_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_M = '%s-%d'%(name_mfn_y_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)

            #* Coordinates
            x = face_nodes[i_node].coordinates[0]

            #* Constraint equations
            '''
            - u1(My) = 0.5*( u1(RP_12) * length_y )
            - u2(My) = 0.5*( u1(RP_12) * x(My)    )
            '''
            myMdl.Equation(name='%s_MY-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 1), (-0.5*length_y, name_rp12, 1)))
            
            myMdl.Equation(name='%s_MY-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 2), (-0.5*x, name_rp12, 1)))

        #* Y-direction: slave face
        
        face_nodes = aa.sets[name_sfn_y_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_S = '%s-%d'%(name_sfn_y_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_S)

            #* Coordinates
            x = face_nodes[i_node].coordinates[0]

            #* Constraint equations
            '''
            - u1(Sy) = 0
            - u2(Sy) = 0.5*( u1(RP_12) * x(Sy)    )
            '''
            myMdl.Equation(name='%s_SY-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 1), (0.0, name_rp12, 1)))
            
            myMdl.Equation(name='%s_SY-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_S, 2), (-0.5*x, name_rp12, 1)))

    @staticmethod
    def create_constraints_epsilon11(myMdl, name_eqn='InPlane11',
                name_mfn_x_set='MFn-X',
                length_x=1.0, name_rp11='RP_11'):
        '''
        Create constraint equations for epsilon11 = u1(RP_11).

        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        name_eqn: str
            name of the equation constraint
            
        name_mfn_x_set: str
            name of the set of master face nodes for x directions
            
        length_x: float
            the length of the model in x direction
            
        name_rp11: str
            name of the reference points for the epsilon11 component

        Notes
        ----------------
        Assuming small deformation, one corner of the RVE is at the origin,
        and the faces goes through the origin are called the slave faces,
        the boundary conditions are:
        
        - u1(Mx) = u1(RP_11) * length_x
        - u2(Mx) = 0
        
        - u1(Sx) = 0
        - u2(Sx) = 0

        '''
        aa = myMdl.rootAssembly
        
        #* X-direction: master face
        
        face_nodes = aa.sets[name_mfn_x_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_M = '%s-%d'%(name_mfn_x_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)

            #* Constraint equations
            '''
            - u1(Mx) = u1(RP_11) * length_x
            - u2(Mx) = 0
            '''
            myMdl.Equation(name='%s_MX-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 1), (-length_x, name_rp11, 1)))
            
            myMdl.Equation(name='%s_MX-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 2), (0.0, name_rp11, 1)))
    
        #* X-direction: slave face
        #* EncastreBC

    @staticmethod
    def create_constraints_epsilon22(myMdl, name_eqn='InPlane22',
                name_mfn_y_set='MFn-Y',
                length_y=1.0, name_rp22='RP_22'):
        '''
        Create constraint equations for epsilon22 = u1(RP_22).

        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        name_eqn: str
            name of the equation constraint
            
        name_mfn_y_set: str
            name of the set of master face nodes for y directions
            
        length_y: float
            the length of the model in y direction
            
        name_rp22: str
            name of the reference point for the epsilon22 component

        Notes
        ----------------
        Assuming small deformation, one corner of the RVE is at the origin,
        and the faces goes through the origin are called the slave faces,
        the boundary conditions are:

        - u1(My) = 0
        - u2(My) = u2(RP_22) * length_y

        - u1(Sy) = 0
        - u2(Sy) = 0
        '''
        aa = myMdl.rootAssembly
        
        #* Y-direction: master face
        
        face_nodes = aa.sets[name_mfn_y_set].nodes

        for i_node in range(len(face_nodes)):

            #* Create set for node
            name_M = '%s-%d'%(name_mfn_y_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((face_nodes[i_node],)), name=name_M)

            #* Constraint equations
            '''
            - u1(My) = 0
            - u2(My) = u2(RP_22) * length_y
            '''
            myMdl.Equation(name='%s_MY-%d-1'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 1), (0.0, name_rp22, 1)))
            
            myMdl.Equation(name='%s_MY-%d-2'%(name_eqn, i_node), 
                    terms=((1.0, name_M, 2), (-length_y, name_rp22, 1)))

        #* Y-direction: slave face
        #* EncastreBC
