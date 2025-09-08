'''
Classes for Abaqus modeling procedures:

-   `PeriodicBC`: class of functions to setup periodic boundary conditions (PBCs) in `Model`;

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


class PeriodicBC(NodeOperation):
    '''
    A class that contains functions to setup periodic boundary conditions (PBCs) in `Model` objects. 
    
    Need to define PBCs after `Mesh`, `Assembly` & `Step`.
    
    PBCs are applied by constraining the relative displacement of some nodes. 
    The constraint equations are introduced to Abaqus by `Constraints -> Equation`.
    
    PBCs are applied between opposite pairs of faces. 
    In each pair of faces, one of them is defined as the master face, the rule to follow is: 
    
        Once a node (that belongs to a master face) is in one equation, this node cannot appear 
        in any other equation.
    
    If we don't follow this rule, we will end up with more equations than unknowns (linear system of equations).
    An attribute `forbiddenNodes` is used to record the nodes that have been constrained.

    Note
    --------------
    If you experience issues with the application of PBCs, this might be because the node pairing is not right. 
    This is a consequence of the precision of the floating point numbers of the coordinates. 
    To solve it, we will round the coordinates of the nodes using a few decimal places (3 or 4 is enough) 
    using the numpy function round: np.round(number, decimal_places).
    
    The equation constraints can be defined by `Constraints -> Equation`, or Abaqus keyword: `*EQUATION` in the input file.
    The keywords can be added to the model using the Keywords Editor.
    But we prefer `Constraints -> Equation` because it's easier to understand and modify in Abaqus CAE.
    
    References
    -------------------    
    [1] Linear constraint equations
    
        A linear multi-point constraint requires that a linear combination of nodal variables is equal to zero.

        https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/usb/default.htm?startat=pt08ch28s02aus104.html#usb-cni-pequation
    
    [2] Defining equation constraints
    
        https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/usi/default.htm?startat=pt03ch15s14hlb07.html
    
        https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/usi/default.htm?startat=pt03ch15s05.html
    
    [3] Abaqus keyword: `*EQUATION`
    
        https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/key/default.htm?startat=ch05abk27.html#usb-kws-mequation
    
    '''
    def __init__(self):
        pass

    @staticmethod
    def exclude_forbidden_nodes_pbc(myMdl, name_instance, master_face_nodes, slave_face_nodes, name_forbidden_sets=[], label_forbidden_nodes=[]):
        '''
        Exclude forbidden nodes from the master face, and exclude the corresponding node from the slave face.
        
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
    def create_node_sets(myMdl, name_instance, name_master_face_set, name_slave_face_set, coords_sorting, name_forbidden_sets=[], label_forbidden_nodes=[]):
        '''
        Create node sets on the master/slave faces for the periodic boundary condition.
        
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
        
        Return
        ----------------
        name_mfn, name_sfn: str
            name of the node sets on the master/slave faces
            
        Examples
        ---------------
        >>> name_mfn, name_sfn = PeriodicBC.create_node_sets( \\
        >>>     myMdl, name_instance, name_master_face_set, name_slave_face_set, 
        >>>     coords_sorting, name_forbidden_sets=[], label_forbidden_nodes=[])
        '''
        #* Sort nodes by coordinates to find node pairs automatically
        master_face_nodes = PeriodicBC.get_nodes_from_face(myMdl, name_master_face_set, coords_sorting, name_instance)
        slave_face_nodes  = PeriodicBC.get_nodes_from_face(myMdl, name_slave_face_set,  coords_sorting, name_instance)
        
        num_node_master = len(master_face_nodes)
        num_node_slave  = len(slave_face_nodes)
        
        if num_node_master != num_node_slave:
            print('>>> --------------------')
            print('[Error]: Periodic boundary conditions (instance: %s)'%(name_instance))
            print('         The number of nodes in the master & slave faces does not match')
            print('         Master face %s (%d); slave face %s (%d)'%(name_master_face_set, num_node_master, name_slave_face_set, num_node_slave))
            raise Exception()
        
        #* Exclude forbidden nodes
        labels_master, labels_slave, label_forbidden_nodes = \
            PeriodicBC.exclude_forbidden_nodes_pbc(myMdl, name_instance, master_face_nodes, slave_face_nodes,
                                                name_forbidden_sets, label_forbidden_nodes)

        name_mfn = 'MFNode-%s-%s'%(name_instance, name_master_face_set)
        name_sfn = 'SFNode-%s-%s'%(name_instance, name_slave_face_set)
        
        #* Create node sets
        myMdl.rootAssembly.SetFromNodeLabels(name=name_mfn, nodeLabels=((name_instance, labels_master),), unsorted=True)
        myMdl.rootAssembly.SetFromNodeLabels(name=name_sfn, nodeLabels=((name_instance, labels_slave ),), unsorted=True)
        
        #* Print information
        print('>>> --------------------')
        print('[Periodic boundary conditions] Instance: %s'%(name_instance))
        print('    Node sets: master face [%s] <=> slave face [%s];'%(name_master_face_set, name_slave_face_set))
        print('    nNodes= %d after excluding forbidden nodes from %d nodes in each face'%(len(labels_master), num_node_master))
        print('>>>')
        
        return name_mfn, name_sfn

    @staticmethod
    def create_constraints_strain_vector(myMdl, name_eqn, name_mfn_set, name_sfn_set, name_mn1, name_mn2,
                                            name_mn3=None, neutral_axis_x=0.0, neutral_axis_y=0.0, 
                                            rotation_axis_x=0.0, rotation_axis_y=0.0):
        '''
        Create constraint equations for specified strain vector in the form of PBC, 
        i.e., [epsilon_z, kappa_z, kappa_x, kappa_y] = [u3_MN1, u1_MN2, u2_MN2, u3_MN2]. 

        Parameters
        ---------------
        myMdl: Abaqus Model object
            model object
            
        name_eqn: str
            name of the equation constraint
            
        name_mfn_set: str
            name of the set of master face nodes
            
        name_sfn_set: str
            name of the set of slave face nodes
            
        name_mn1, name_mn2: str
            name of the set of the master node 1 & 2

        name_mn3: None, or str
            If the neutral axis of bending is to be solved, provide name of the set of the master node 3.

        neutral_axis_x, neutral_axis_y: float
            the coordinates of the neutral axis of bending.
            
        rotation_axis_x, rotation_axis_y: float
            The rotation axis is parallel to z-axis and foes through (rotation_axis_x, rotation_axis_y)
            
        Notes
        ----------------       
        Assuming small deformation, the boundary conditions are implemented using master node `MN1` and `MN2`:
        
        - `u1_M` - `u1_S` = `u1_MN1` - (`y`-`rotation_axis_y`)*`L`*`u1_MN2`
        - `u2_M` - `u2_S` = `u2_MN1` + (`x`-`rotation_axis_x`)*`L`*`u1_MN2`
        - `u3_M` - `u3_S` = `L`*`u3_MN1` + (`y`-`neutral_axis_y`)*`L`*`u2_MN2` - (`x`-`neutral_axis_x`)*`L`*`u3_MN2`

        It is assumed that the master and slave faces are approximately in the x-y plane.
        
        - The tension in    z-axis is applied by giving `u3_MN1` a displacement corresponding to a desired axial strain.
        - The torsion about z-axis is applied by giving `u1_MN2` a displacement corresponding to a desired twist ratio,
        - The bending about x-axis is applied by giving `u2_MN2` a displacement corresponding to a desired curvature.
        - The bending about y-axis is applied by giving `u3_MN2` a displacement corresponding to a desired curvature.
        
        - The relative rigid displacement between the master and slave faces in the x-y plane is not constrained, i.e., `u1_MN1`, `u2_MN1` are left free.
        
        The variables in the functions are:
        
        - `u` is the displacement of nodes on the master and slave faces;
        - `L` is the length of model in z-axis direction;
        - `x` is the distance from the rotation axis (, also the neutral axis of bending about y-axis);
        - `y` is the distance from the rotation axis (, also the neutral axis of bending about x-axis);
        
        If the neutral axis of bending is to be solved, the 3rd constraint equation becomes:
        
        - `u3_M` - `u3_S` = `L`*`u3_MN1` + `y`*`L`*`u2_MN2` - `x`*`L`*`u3_MN2` + `u3_MN3`
        
        where,
        
        - `u3_MN3` = `neutral_axis_y`*`L`*`u2_MN2` + `neutral_axis_x`*`L`*`u3_MN2`
        
        References
        ------------------
        R. Hasa, S.T. Pinho, Failure mechanisms of biological crossed-lamellar micro-structures applied to 
        synthetic high-performance fibre-reinforced composites, Journal of the Mechanics and Physics of Solids, 2019
        '''
        aa = myMdl.rootAssembly
        
        master_face_nodes = aa.sets[name_mfn_set].nodes
        slave_face_nodes  = aa.sets[name_sfn_set].nodes
        
        n_node = len(master_face_nodes)
        
        for i_node in range(n_node):

            #* Create sets for each pair of nodes
            name_M = '%s-%d'%(name_mfn_set, i_node)
            name_S = '%s-%d'%(name_sfn_set, i_node)
            aa.Set(nodes=mesh.MeshNodeArray((master_face_nodes[i_node],)), name=name_M)
            aa.Set(nodes=mesh.MeshNodeArray((slave_face_nodes [i_node],)), name=name_S)

            #* Coefficients
            x = master_face_nodes[i_node].coordinates[0]
            y = master_face_nodes[i_node].coordinates[1]

            #* Constraint equations
            myMdl.Equation(name='%s-%d-x'%(name_eqn, i_node), terms=((1.0, name_M, 1), (-1.0, name_S, 1), (-1.0, name_mn1, 1),
                                                    ( (y-rotation_axis_y), name_mn2, 1)))
            
            myMdl.Equation(name='%s-%d-y'%(name_eqn, i_node), terms=((1.0, name_M, 2), (-1.0, name_S, 2), (-1.0, name_mn1, 2),
                                                    (-(x-rotation_axis_x), name_mn2, 1)))
            
            if name_mn3 is None:
            
                myMdl.Equation(name='%s-%d-z'%(name_eqn, i_node), terms=((1.0, name_M, 3), (-1.0, name_S, 3), (-1, name_mn1, 3), 
                                                    (-(y-neutral_axis_y), name_mn2, 2), ( (x-neutral_axis_x), name_mn2, 3)))

            else:
                
                myMdl.Equation(name='%s-%d-z'%(name_eqn, i_node), terms=((1.0, name_M, 3), (-1.0, name_S, 3), (-1, name_mn1, 3), 
                                                    (-y, name_mn2, 2), ( x, name_mn2, 3),
                                                    ( 1.0, name_mn3, 3)))
    
    @staticmethod
    def calculate_master_node_displacement_BC(strain_vector, model_length_z):
        '''
        Calculate the displacement `u` of master nodes for boundary condition defination.

        Parameters
        ----------------------
        strain_vector: ndarray [4]
            the strain vector for loading
            
        model_length_z: float
            the length of model `L` in z-axis direction
        
        Returns
        ----------------------
        u3_MN1: float
            the displacement corresponding to a desired axial strain for tension in z-axis
        
        u1_MN2: float
            the displacement corresponding to a desired twist ratio for torsion about z-axis
            
        u2_MN2: float
            the displacement corresponding to a desired curvature for bending about x-axis
        
        u3_MN2: float
            the displacement corresponding to a desired curvature for bending about y-axis
        '''
        u3_MN1 = strain_vector[0]*model_length_z
        u1_MN2 = strain_vector[1]*model_length_z
        u2_MN2 = strain_vector[2]*model_length_z
        u3_MN2 = strain_vector[3]*model_length_z

        return u3_MN1, u1_MN2, u2_MN2, u3_MN2
           
    @staticmethod
    def calculate_bending_x_neutral_axis(u3_MN3, u2_MN2, model_length_z):
        '''
        Calculate the neutral axis of bending about x-axis with the displacement of master nodes.
        The result is only valid when there is no rigid motion for the pure bending case.
        
        Parameters
        ---------------
        u3_MN3: float
            the displacement of master node 3
        
        u2_MN2: float
            the displacement of master node 2
            
        model_length_z: float
            the length of model `L` in z-axis direction.
            
        Returns
        ---------------
        yB: float
            the coordinate of the neutral axis
        '''
        yB = u3_MN3 / u2_MN2 / model_length_z
        return yB

    @staticmethod
    def calculate_bending_y_neutral_axis(u3_MN3, u3_MN2, model_length_z):
        '''
        Calculate the neutral axis of bending about y-axis with the displacement of master nodes.
        The result is only valid when there is no rigid motion for the pure bending case.

        Parameters
        ---------------
        u3_MN3: float
            the displacement of master node 3
        
        u3_MN2: float
            the displacement of master node 2
            
        model_length_z: float
            the length of model `L` in z-axis direction.
            
        Returns
        ---------------
        xB: float
            the coordinate of the neutral axis
        '''
        xB = u3_MN3 / u3_MN2 / model_length_z
        return xB

