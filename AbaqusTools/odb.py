'''
Classes for Abaqus modeling procedures:

-   `Part`: class for creating an Abaqus part,
    including procedures in *Sketch*, *Part*, *Property* and *Mesh* modules;
-   `Model`: class for creating an Abaqus model,
    including model initialization, building `Part` s, and procedures in
    *Assembly*, *Step*, *Interaction*, *Load* and *Job* modules;
-   `PeriodicBC`: class of functions to setup periodic boundary conditions (PBCs) in `Model`;
'''
import numpy as np

from AbaqusTools import IS_ABAQUS

if IS_ABAQUS:

    from abaqus import *
    from abaqusConstants import *
    from symbolicConstants import *
    from caeModules import *
    import mesh
    import odbAccess


class OdbOperation(object):
    '''
    The operations to the output database of a specified job.
    
    Data is stored and calculated in the integration points of elements.
    The value on nodal is interpolated from the integration points.
    
    Reference
    ----------------
    https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo04.html
    
    '''
    def __init__(self, name_job):
        
        self.name_job = name_job+'.odb'
        
        #* Mapping between the fieldValue index and the label in its instance
        self.mapping_node_label2index = None
        self.mapping_node_index2label = None
        
        self.mapping_elem_label2index = None
        self.mapping_elem_index2label = None
        
        #* Loads odb file and create a new Odb object
        session.openOdb(name=self.name_job, path=self.name_job)
    
    @property
    def odb(self):
        '''
        Output database of job
        
        Odb object:
        
        https://docs.software.vt.edu/abaqusv2022/English/?show=SIMACAEKERRefMap/simaker-c-odbpyc.htm
        '''
        return session.odbs[self.name_job]
    
    def get_instance(self, name_instance):
        '''
        Get an instance of the output database.
        
        OdbInstance object:
        
            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html
        
        Parameters
        ---------------
        name_instance: str
            name of an instance in CAPITAL letters, e.g., 'LOCAL_SPAR_FRONT'.
            
        Returns
        ---------------
        instance: OdbInstance
            a part instance is the usage of a part within an assembly, which has the following members:
            - name: a string specifying the instance name.;
            - nodes: an OdbMeshNodeArray object;
            - elements: an OdbMeshElementArray object;
            - nodeSets: a repository of OdbSet objects specifying node sets;
            - elementSets: a repository of OdbSet objects specifying element sets;
            - surfaces: a repository of OdbSet objects specifying surfaces.
        '''
        return session.odbs[self.name_job].rootAssembly.instances[name_instance]
    
    def get_nodes(self, name_instance, name_set=None, name_surface=None):
        '''
        Get nodes in a set of an instance.
        
        FieldValue object:
        
            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html
        
        OdbMeshNode object:
        
            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html
        
        Parameters
        ---------------
        name_instance: str
            name of an instance in CAPITAL letters, e.g., 'SPAR_FRONT'.
            
        name_set: None, or str
            name of a set in CAPITAL letters, e.g., 'LOCAL-FACE-WEB-U'.
            Default is None, which outputs all nodes in the instance.
            
        name_surface: None, or str
            name of a surface in CAPITAL letters, e.g., 'BOLT-1'.
            Default is None, which outputs all nodes in the instance.
        
        Returns
        ---------------
        nodes: OdbMeshNodeArray
            a list of OdbMeshNode object, which has the following members:
            - label: an Int specifying the node label;
            - coordinates: a tuple of Floats specifying the nodal coordinates in the global Cartesian coordinate system.
        '''
        instance = self.get_instance(name_instance)
        
        if name_set is None and name_surface is None:
            return instance.nodes
        
        elif name_set is not None and name_surface is not None:
            print('>>> --------------------')
            print('Error [OdbOperation.get_nodes]: ')
            print('    Do not input [name_set] and [name_surface] at the same time: [%s], [%s]', name_set, name_surface)
            raise Exception
        
        elif name_set is not None:
            return instance.nodeSets[name_set].nodes
        
        elif name_surface is not None:
            return instance.surfaces[name_surface].nodes
        
        else:
            raise Exception
    
    def get_elements(self, name_instance, name_set=None, name_surface=None):
        '''
        Get elements in a set of an instance.
        
        FieldValue object:
        
            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html
        
        OdbMeshElement object:
        
            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html
        
        Parameters
        ---------------
        name_instance: str
            name of an instance in CAPITAL letters, e.g., 'LOCAL_SPAR_FRONT'.
            
        name_set: None, or str
            name of a set in CAPITAL letters, e.g., 'FACE_BOLT_0'.
            Default is None, which outputs all nodes in the instance.
        
        name_surface: None, or str
            name of a surface in CAPITAL letters, e.g., 'BOLT-1'.
            Default is None, which outputs all nodes in the instance.
        
        Returns
        ---------------
        elements: OdbMeshElementArray
            a list of OdbMeshElement object, which has the following members:
            - label: an Int specifying the element label;
            - type: a string specifying the element type;
            - sectionCategory: a SectionCategory object specifying the element section properties;
            - connectivity: a tuple of Ints specifying the element connectivity;
            - instanceNames: a tuple of Strings specifying the instance names for nodes in the element connectivity.

        '''
        instance = self.get_instance(name_instance)
        
        if name_set is None and name_surface is None:
            return instance.elements
        
        elif name_set is not None and name_surface is not None:
            print('>>> --------------------')
            print('Error [OdbOperation.get_elements]: ')
            print('    Do not input [name_set] and [name_surface] at the same time: [%s], [%s]', name_set, name_surface)
            raise Exception
        
        elif name_set is not None:
            return instance.elementSets[name_set].elements
        
        elif name_surface is not None:
            return instance.surfaces[name_surface].elements
        
        else:
            raise Exception

    def get_fieldOutput(self, step='Loading', frame=-1, variable='U'):
        '''
        The FieldOutput object of the last frame of the 'Loading' step.

        FieldOutput object:
        
            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo04.html
        
            https://docs.software.vt.edu/abaqusv2022/English/?show=SIMACAEKERRefMap/simaker-c-fieldoutputpyc.htm
        
        FieldLocation object:
        
            The FieldLocation object specifies locations for which data are available in the field.
            
            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo03.html
        
        Parameters
        ----------------
        step: str
            name of the step, e.g., 'Loading'
        
        frame: int
            index of the frame, default -1 means the last frame
            
        variable: str
            name of the output variable, e.g., 'S', 'E', 'U', 'RF', etc.
        
        Return
        ---------
        fieldOutput: FieldOutput object
            It contains the following members:
            - componentLabels: a tuple of strings, e.g., ('RF1', 'RF2', 'RF3')
            - description: a string, e.g., 'Reaction force'
            - locations: a FieldLocationArray object
            - values: a FieldValueArray object
            
        position: str
            The text of a SymbolicConstant specifying the position of the output in the element. Possible values are:
            - NODAL: specifying the values calculated at the nodes.
            - INTEGRATION_POINT: specifying the values calculated at the integration points.
            - ELEMENT_NODAL: specifying the values obtained by extrapolating results calculated at the integration points.
            - CENTROID: specifying the value at the centroid obtained by extrapolating results calculated at the integration points.
            
        data_type: str
            The text of a SymbolicConstant specifying output type.
            Possible values are SCALAR, VECTOR, TENSOR_3D_FULL, TENSOR_3D_PLANAR, TENSOR_3D_SURFACE, TENSOR_2D_PLANAR, and TENSOR_2D_SURFACE.
        '''
        fieldOutput = session.odbs[self.name_job].steps[step].frames[frame].fieldOutputs[variable]
        
        position = fieldOutput.locations[0].position.getText()
        
        data_type = fieldOutput.type.getText()
        
        return fieldOutput, position, data_type
    
    def get_num_frames(self, step='Loading'):
        '''
        Get the number of frames in a step.
        
        Parameters
        ----------------
        step: str
            name of the step, e.g., 'Loading'
        '''
        return len(session.odbs[self.name_job].steps[step].frames)
    
    @staticmethod
    def get_fieldValue(fieldOutput, index):
        '''
        Get a FieldValue object from a FieldOutput object,
        which represents the field data at a point.
        
        FieldValue object:
        
            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html

        Parameters
        ----------------
        fieldOutput: FieldOutput object
            it contains field data for a specific output variable.
        
        Return
        ---------
        fieldValue: FieldValue object
            It contains the following members:
            - position: a SymbolicConstant specifying the position of the output in the element. Possible values are: NODAL, INTEGRATION_POINT, ELEMENT_NODAL, ELEMENT_FACE, CENTROID.
            - elementLabel: an Int specifying the element label of the element containing the location. elementLabel is available only if position=INTEGRATION_POINT, CENTROID, ELEMENT_NODAL, or ELEMENT_FACE.
            - nodeLabel: an Int specifying the node label of the node containing the location. nodelabel is available only if position=ELEMENT_NODAL or NODAL.
            - integrationPoint: An Int specifying the integration point in the element. integrationPoint is available only if position=INTEGRATION_POINT.
            - data, dataDouble: a tuple of Floats specifying data
        
        type_str: str
            the text for a SymbolicConstant specifying the output type. 
            Possible values are SCALAR, VECTOR, TENSOR_3D_FULL, TENSOR_3D_PLANAR, TENSOR_3D_SURFACE, TENSOR_2D_PLANAR, and TENSOR_2D_SURFACE.
        
        data: ndarray
            a 1-d array of values
        '''
        fieldValue = fieldOutput.values[index]
        type_str = fieldValue.type.getText()
        data = fieldValue.data
        
        if type_str == 'SCALAR':
            data = np.array([data])
            
        elif type_str == 'VECTOR':
            data = np.array(data)
            
        elif type_str == 'TENSOR_3D_FULL':
            data = np.array(data)
            
        else:
            raise NotImplementedError
            
        return fieldValue, type_str, data
    
    #* Get node/element's label/index
    def create_node_index_mapping(self, variable='U'):
        '''
        Create the mapping between the node label in its instance and the node index in the fieldOutput.
        
        Parameters
        ----------------
        variable: str
            name of the output variable, e.g., 'U', etc.
        
        Attributes
        ----------------
        mapping_node_label2index: dict of a dict
        
            Its key is the name of instance in CAPITAL letters, e.g., 'SPAR_FRONT'.
            
            Its value is a dictionary. The secondary dictionary's key is the node label in its instance,
            the value is its index in the fieldOutput.
            
            Therefore, index_fieldOutput = mapping_node_label2index[name_instance][node_label].
        
        mapping_node_index2label: list of tuple
            A list of (name_instance, node_label)
        '''
        self.mapping_node_label2index = {}
        self.mapping_node_index2label = []
        
        fieldOutput, position, _ = self.get_fieldOutput(variable=variable)
        
        if not position == 'NODAL':
            print('>>> --------------------')
            print('Error [OdbOperation.create_node_index_mapping]: ')
            print('    Variable [%s] is stored in [%s] instead of NODAL'%(variable, position))
            raise Exception
        
        # fieldOutput.values is a list of FieldValue objects
        n_total = len(fieldOutput.values)
        
        for i_fieldOutput in range(n_total):

            # FieldValue object
            # https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html
            
            fieldValue = fieldOutput.values[i_fieldOutput]
            
            # Get the name of instance
            if fieldValue.instance is not None:
                name_instance = fieldValue.instance.name
            else:
                name_instance = 'ASSEMBLY'
                
            if name_instance not in self.mapping_node_label2index.keys():
                self.mapping_node_label2index[name_instance] = {}
            
            # Get the label of node in its instance, it starts from 1
            node_label = fieldValue.nodeLabel
            self.mapping_node_label2index[name_instance][node_label] = i_fieldOutput
            
            self.mapping_node_index2label.append((name_instance, node_label))
        
        print('>>> [OdbOperation]: %s'%(self.name_job))
        print('         The mapping between node label and its index in fieldOutput is created.')
        print('>>>')
    
    def create_element_index_mapping(self, variable='S'):
        '''
        Create the mapping between the element label in its instance and the element index in the fieldOutput.
        
        Parameters
        ----------------
        variable: str
            name of the output variable, e.g., 'S', 'E', 'UVARM1', etc.
        
        Attributes
        ----------------
        mapping_elem_label2index: dict of a dict
        
            Its key is the name of instance in CAPITAL letters, e.g., 'LOCAL_SPAR_FRONT'.
            
            Its value is a dictionary. The secondary dictionary's key is the element label in its instance,
            the value is its index in the fieldOutput.
            
            Therefore, index_fieldOutput = mapping_elem_label2index[name_instance][element_label].
        
        mapping_elem_index2label: list of tuple
            A list of (name_instance, element_label)
        '''
        self.mapping_elem_label2index = {}
        self.mapping_elem_index2label = []
        
        fieldOutput, position, _ = self.get_fieldOutput(variable=variable)
        
        if not position == 'INTEGRATION_POINT':
            print('>>> --------------------')
            print('Error [OdbOperation.create_element_index_mapping]: ')
            print('    Variable [%s] is stored in [%s] instead of INTEGRATION_POINT'%(variable, position))
            raise Exception
        
        # fieldOutput.values is a list of FieldValue objects
        n_total = len(fieldOutput.values)
        
        for i_fieldOutput in range(n_total):

            # FieldValue object
            # https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html
            
            fieldValue = fieldOutput.values[i_fieldOutput]
            
            # Get the name of instance
            if fieldValue.instance is not None:
                name_instance = fieldValue.instance.name
            else:
                name_instance = 'ASSEMBLY'
                
            if name_instance not in self.mapping_elem_label2index.keys():
                self.mapping_elem_label2index[name_instance] = {}
            
            # Get the label of element in its instance, it starts from 1
            element_label = fieldValue.elementLabel
            self.mapping_elem_label2index[name_instance][element_label] = i_fieldOutput
            
            self.mapping_elem_index2label.append((name_instance, element_label))
        
        print('>>> [OdbOperation]: %s'%(self.name_job))
        print('         The mapping between element label and its index in fieldOutput is created.')
        print('>>>')
    
    def get_node_labels_and_indices(self, name_instance, name_set=None, name_surface=None):
        '''
        Get the the labels and indices of nodes in a set of an instance.
        
        https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html
        
        Parameters
        ---------------
        name_instance: str
            name of an instance in CAPITAL letters, e.g., 'SPAR_FRONT'.
            
        name_set: None, or str
            name of a set in CAPITAL letters, e.g., 'LOCAL-FACE-WEB-U'.
            Default is None, which outputs all nodes in the instance.
        
        name_surface: None, or str
            name of a surface in CAPITAL letters, e.g., 'BOLT-1'.
            Default is None, which outputs all nodes in the instance.
        
        Returns
        ---------------
        node_labels: list of int
            labels of nodes in the instance, starts from 1.
            It equals the node index in the OdbMeshNodeArray `nodes` plus 1.
        
        indices_fieldOutput: list of int
            indices of nodes in the fieldOutput
        '''
        nodes = self.get_nodes(name_instance, name_set, name_surface)
        
        if self.mapping_node_label2index is None:
            self.create_node_index_mapping()
        
        node_labels = []
        indices_fieldOutput = []
        
        for i in range(len(nodes)):
            
            node_labels.append(nodes[i].label)
            
            indices_fieldOutput.append(
                self.mapping_node_label2index[name_instance][nodes[i].label])

        return node_labels, indices_fieldOutput
    
    def get_element_labels_and_indices(self, name_instance, name_set=None, name_surface=None):
        '''
        Get the the labels and indices of elements in a set of an instance.
        
        https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo05.html
        
        Parameters
        ---------------
        name_instance: str
            name of an instance in CAPITAL letters, e.g., 'LOCAL_SPAR_FRONT'.
            
        name_set: None, or str
            name of a set in CAPITAL letters, e.g., 'FACE_BOLT_0'.
            Default is None, which outputs all elements in the instance.
        
        name_surface: None, or str
            name of a surface in CAPITAL letters, e.g., 'BOLT-1'.
            Default is None, which outputs all nodes in the instance.
        
        Returns
        ---------------
        element_labels: list of int
            labels of elements in the instance, starts from 1.
            It equals the element index in the OdbMeshElementArray `elements` plus 1.
        
        indices_fieldOutput: list of int
            indices of elements in the fieldOutput
        '''
        elements = self.get_elements(name_instance, name_set, name_surface)
        
        if self.mapping_elem_label2index is None:
            self.create_element_index_mapping()
        
        element_labels = []
        indices_fieldOutput = []
        
        for i in range(len(elements)):
            
            element_labels.append(elements[i].label)
            
            indices_fieldOutput.append(
                self.mapping_elem_label2index[name_instance][elements[i].label])

        return element_labels, indices_fieldOutput
    
    #* Probe value of node/element
    def probe_node_values(self, step='Loading', frame=-1, variable='U', component=None, index_fieldOutput=0):
        '''
        Probe values of a node or nodes.
        
        Parameters
        ----------------
        step: str
            name of the step, e.g., 'Loading'
        
        frame: int
            index of the frame, default -1 means the last frame
            
        variable: str
            name of the output variable, e.g., 'S', 'U', etc.
        
        component: None, or str
            name of the component, e.g., 'U1', 'S11', etc.
            Default None means output all components.
        
        index_fieldOutput: int or list[int]
            index of the node in the fieldOutputs, which contains all nodes from all instances.
        
        Returns
        ---------------
        values: ndarray [n_comp] or [n_node, n_comp]
            an array of values
        '''
        fieldOutput, position, _ = self.get_fieldOutput(step, frame, variable)

        if not position == 'NODAL':
            
            print('>>> --------------------')
            print('Error [probe_node_values]: the variable is not stored in nodes')
            print('    Step: [%s]; Frame: [%d]'%(step, frame))
            print('    The location of field data for [%s] is [%s]'%(variable, position))
            raise Exception()
        
        if isinstance(index_fieldOutput, int):

            data = fieldOutput.values[index_fieldOutput].data

            if component is None:
                
                values = data
            
            else:
            
                index_comp = fieldOutput.componentLabels.index(component)
                values = data[index_comp:index_comp+1]
                
        else:
            
            values = []
            
            for i_node in range(len(index_fieldOutput)):
                
                data = fieldOutput.values[index_fieldOutput[i_node]].data
            
                if component is None:
                    
                    values.append(data)
                
                else:
                
                    index_comp = fieldOutput.componentLabels.index(component)
                    values.append(data[index_comp:index_comp+1])
            
            values = np.array(values)
            
        return values
    
    def probe_element_values(self, step='Loading', frame=-1, variable='U', component=None, index_fieldOutput=0):
        '''
        Probe values of a element or elements. The value is stored in integration point(s).
        
        Parameters
        ----------------
        step: str
            name of the step, e.g., 'Loading'
        
        frame: int
            index of the frame, default -1 means the last frame
            
        variable: str
            name of the output variable, e.g., 'S', 'U', etc.
        
        component: None, or str
            name of the component, e.g., 'U1', 'S11', etc.
            Default None means output all components.
        
        index_fieldOutput: int or list[int]
            index of the element in the fieldOutputs, which contains all elements from all instances.
        
        Returns
        ---------------
        values: ndarray [n_comp] or [n_node, n_comp]
            an array of values
        '''
        fieldOutput, position, _ = self.get_fieldOutput(step, frame, variable)

        if not position == 'INTEGRATION_POINT':
            
            print('>>> --------------------')
            print('Error [probe_element_values]: the variable is not stored in elements')
            print('    Step: [%s]; Frame: [%d]'%(step, frame))
            print('    The location of field data for [%s] is [%s]'%(variable, position))
            raise Exception()
        
        if isinstance(index_fieldOutput, int):

            data = fieldOutput.values[index_fieldOutput].data

            if component is None:
                
                values = data
            
            else:
            
                index_comp = fieldOutput.componentLabels.index(component)
                values = data[index_comp:index_comp+1]
                
        else:
            
            values = []
            
            for i_node in range(len(index_fieldOutput)):
                
                data = fieldOutput.values[index_fieldOutput[i_node]].data
            
                if component is None:
                    
                    values.append(data)
                
                else:
                
                    index_comp = fieldOutput.componentLabels.index(component)
                    values.append(data[index_comp:index_comp+1])
            
            values = np.array(values)
            
        return values
    
    #* Probe coordinate of node/element
    def probe_node_coordinate(self, name_instance='ASSEMBLY', node_label=1, index_fieldOutput=None):
        '''
        Get the coordinate of a node or nodes. 
        Need to provide either the node label in its instance, or its index in the fieldOutput.
        
        OdbMeshNodeArray object:

            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo18.html
        
        Parameters
        --------------
        name_instance: str
            name of an instance in CAPITAL letters, e.g., 'ASSEMBLY'.
        
        node_label: int or list[int]
            label of node in the instance, it starts from 1.
            It equals the node index in the OdbMeshNodeArray `nodes` plus 1.
        
        index_fieldOutput: int or list[int]
            index of the node in the fieldOutputs, which contains all nodes from all instances.
        
        Returns
        ---------------
        coordinates: ndarray [n_dim] or [n_node, n_dim]
            the nodal coordinates in the global Cartesian coordinate system
        '''
        mesh_nodes = session.odbs[self.name_job].rootAssembly.instances[name_instance].nodes
        
        #* Convert index to label
        if index_fieldOutput is not None:
            
            if self.mapping_node_index2label is None:
                self.create_node_index_mapping()
            
            if isinstance(index_fieldOutput, int):

                name_instance, node_label = self.mapping_node_index2label[index_fieldOutput]
                    
            else:
                
                node_label = []
                for i_node in range(len(index_fieldOutput)):
                    
                    _, label = self.mapping_node_index2label[index_fieldOutput[i_node]]
                    node_label.append(label)
        
        #* Get coordinates
        if isinstance(node_label, int):

            coordinates = np.array(mesh_nodes[node_label-1].coordinates)
            
        else:
            
            coordinates = []
            for i_node in range(len(node_label)):
                coordinates.append(mesh_nodes[node_label[i_node]-1].coordinates)
            
            coordinates = np.array(coordinates)

        return coordinates

    def probe_element_center_coordinate(self, name_instance='ASSEMBLY', element_label=1, index_fieldOutput=None):
        '''
        Get the coordinate of the center of an element. 
        Need to provide either the element label in its instance, or its index in the fieldOutput.
        
        Parameters
        --------------
        name_instance: str
            name of an instance in CAPITAL letters, e.g., 'ASSEMBLY'.
        
        element_label: int or list[int]
            label of element in the instance, it starts from 1.
            It equals the element index in the OdbMeshElementArray `elements` plus 1.
        
        index_fieldOutput: int or list[int]
            index of the element in the fieldOutputs, which contains all elements from all instances.
        
        Returns
        ---------------
        coordinates: ndarray [n_dim] or [n_node, n_dim]
            the nodal coordinates in the global Cartesian coordinate system
        '''
        #* Get connectivity of element(s)
        connectivity, instanceNames, element_label = self.get_element_connectivity(name_instance, element_label, index_fieldOutput)
        
        n_element = len(element_label)
        if n_element == 1:
            connectivity = [connectivity]
            instanceNames = [instanceNames]
        
        #* Get coordinates
        coordinates = []
        
        for i_elem in range(n_element):

            for i_node in range(len(instanceNames[i_elem])):
                name = instanceNames[i_elem][i_node]
                if not name == name_instance:
                    print('>>> Error [OdbOperation.probe_element_center_coordinate]')
                    print('    Element instance: '+name_instance)
                    print('    Element label:    %d'%(element_label[i_elem]))
                    print('    Node instance:    '+name)
                    print('    Node label:       %d'%(connectivity[i_elem][i_node]))
                    raise Exception()

            XYZs = self.probe_node_coordinate(name_instance, node_label=connectivity[i_elem])
            XYZs = np.mean(XYZs, axis=0)
            
            coordinates.append(XYZs)

        coordinates = np.array(coordinates)
        if n_element == 1:
            coordinates = coordinates[0,:]
        
        return coordinates

    def get_element_connectivity(self, name_instance='ASSEMBLY', element_label=1, index_fieldOutput=None):
        '''
        Get the connectivity of an element, i.e., the label of nodes connected to the element. 
        Need to provide either the element label in its instance, or its index in the fieldOutput.
        
        OdbMeshNodeArray object:

            https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/ker/default.htm?startat=pt01ch31pyo18.html
        
        Parameters
        --------------
        name_instance: str
            name of an instance in CAPITAL letters, e.g., 'ASSEMBLY'.
        
        element_label: int or list[int]
            label of element in the instance, it starts from 1.
            It equals the element index in the OdbMeshElementArray `elements` plus 1.
        
        index_fieldOutput: int or list[int]
            index of the element in the fieldOutputs, which contains all elements from all instances.
        
        Returns
        ---------------
        connectivity: list[int] or list[list[int]]
            the label of nodes connected to the element. 
            
        instanceNames: list[str] or list[list[str]]
            the names of instance that the node belongs to.
            
        element_label: list[int]
            the label of element(s)
        '''
        #* Convert index to label
        if index_fieldOutput is not None:
            
            if self.mapping_elem_index2label is None:
                self.create_element_index_mapping()
            
            if isinstance(index_fieldOutput, int):

                name_instance, element_label = self.mapping_elem_index2label[index_fieldOutput]
                    
            else:
                
                element_label = []
                for i_elem in range(len(index_fieldOutput)):
                    
                    _, label = self.mapping_elem_index2label[index_fieldOutput[i_elem]]
                    element_label.append(label)

        #* Get connectivity
        instance = self.get_instance(name_instance)
        
        if isinstance(element_label, int):
            
            element = instance.elements[element_label-1]
            
            connectivity  = list(element.connectivity)
            instanceNames = list(element.instanceNames)
            element_label = [element_label]
            
        else:
            
            connectivity  = []
            instanceNames = []

            for i_elem in range(len(element_label)):
                
                element = instance.elements[element_label[i_elem]-1]
                connectivity.append(list(element.connectivity))
                instanceNames.append(list(element.instanceNames))
        
        return connectivity, instanceNames, element_label


