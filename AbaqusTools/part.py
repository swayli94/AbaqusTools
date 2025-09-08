'''
Classes for Abaqus modeling procedures:

-   `Part`: class for creating an Abaqus part,
    including procedures in *Sketch*, *Part*, *Property* and *Mesh* modules;

'''
import numpy as np
from scipy.optimize import fsolve
import time

from AbaqusTools import IS_ABAQUS

if IS_ABAQUS:

    from abaqus import *
    from abaqusConstants import *
    from symbolicConstants import *
    from caeModules import *
    import mesh
    import odbAccess


class Part(object):
    '''
    Class for creating an Abaqus Part, 
    including operations in the `Sketch`, `Part`, `Property` and `Mesh` Abaqus Modules.
    
    Parameters
    --------------
    model: Abaqus model object, or None
        Abaqus model. It can be None for testing.
        
    pGeo: dict, or None
        geometry parameters. It can be None for testing.
        
    pMesh: dict, or None
        meshing parameters. It can be None for testing.
    '''
    def __init__(self, model, pGeo=None, pMesh=None):
        
        self.model = model
        self.pGeo = pGeo
        self.pMesh = pMesh
        
        self.raise_exception = True
        self.name_part = 'Part'
        
        self.is_only_geometry = False
    
    #* =============================================
    #* Abaqus procedure
    def build(self):
        '''
        Build an Abaqus part:
        
        - operations in the `Sketch` Abaqus Modules
        
        - operations in the `Part` Abaqus Modules
        
            - Create Part
            - Create Partition
            - Create surfaces and sets
            
        - operations in the `Mesh` Abaqus Modules
        
            - Set seeding
            - Create Mesh
            - Assign Element Type
            
        - operations in the `Property` Abaqus Modules
        
            - Assign Section
            - Create Composite Layup
        '''
        t0 = time.time()
        
        #* Abaqus Module: Sketch
        self.create_sketch()
        
        #* Abaqus Module: Part
        self.create_part()
        self.create_partition()
        self.create_surface()
        self.create_set()
        
        if not self.is_only_geometry:
        
            #* Abaqus Module: Mesh
            self.set_seeding()
            self.create_mesh()
            self.set_element_type()
            
            #* Abaqus Module: Property
            self.set_section_assignment()
            self.set_composite_layups()
        
        t1 = time.time()

        print('>>> --------------------')
        print('    [Part: %s] build time = %.1f (min)'%(self.name_part, (t1-t0)/60.0))
        print('>>>')
        
    def create_sketch(self):
        '''
        Create stand-alone sketches via Abaqus Module: Sketch
        
        - Create Sketch
        '''
    
    def create_part(self):
        '''
        Create part via Abaqus Module: Part
        
        - Create Datum Plane/Axis/CSYS
        - Create Part
        - Create Solid: Extrude
        - Create Shell: Extrude
        - Create Cut: Extrude
        - Create Round and Fillet
        - Create Mirror
        '''

    def create_partition(self):
        '''
        Create partition via Abaqus Module: Part
        
        - Create Partition
        - Partition Face: Sketch
        - Partition Face: Use Datum Plane
        - Partition Cell: Define Cutting Plane
        - Partition Cell: Use Datum Plane
        - Partition Cell: Extend Face
        - Partition Cell: Extrude/Sweep Edges
        '''

    def create_surface(self):
        '''
        Create surfaces
        '''

    def create_set(self):
        '''
        Create sets
        '''

    def set_seeding(self):
        '''
        Set seeding via Abaqus Module: Mesh
        
        - Seed Edges
        - Seed Part
        '''

    def create_mesh(self):
        '''
        Set mesh control and create mesh via Abaqus Module: Mesh
        
        - Assign Mesh Control
        - Assign Stack Direction
        - Mesh Region
        - Mesh Part
        '''
    
    def set_element_type(self):
        '''
        Set element type via Abaqus Module: Mesh
        
        - Assign Element Type
        '''
    
    def set_section_assignment(self):
        '''
        Set section assignment via Abaqus Module: Property
        
        - Assign Section
        '''
    
    def set_composite_layups(self):
        '''
        Set composite layup via Abaqus Module: Property
        
        - Create Composite Layup
        
        References
        ---------------------
        Creating a composite layup
        
            http://130.149.89.49:2080/v6.10/books/usi/default.htm?startat=pt04ch22.html
        
        Creating continuum shell composite layups
        
            http://130.149.89.49:2080/v6.10/books/usi/default.htm?startat=pt03ch12s13s03.html
        '''

    #* =============================================
    #* Abaqus Part functions
    @property
    def part(self, name_part):
        '''
        Abaqus Part object by the name of `name_part` (str).
        '''
        return self.model.parts[name_part]
        
    def create_geometry_set(self, name_set, name_part, findAt_points, 
                            geometry='vertex', getClosest=True, searchTolerance=1E-6):
        '''
        Create a geometry set in Part
        
        Parameters
        ----------------
        name_set: str
            name of the set
            
        name_part: str
            name of the part that the geometry belong to.
            
        findAt_points: tuple, or a list of tuple
        
            point coordinates.
            
            Each face is specified by one point coordinate tuple (x,y,z).
            
            `findAt_points` can be either a tuple of one point, or a list of point tuples.
            
        geometry: str
            'vertex', 'edge', 'face', 'cell'
            
        getClosest: bool
            whether use `getClosest` to find an vertex or edge
            
        searchTolerance: float
            the distance within which the closest object must lie
        '''
        myPrt = self.model.parts[name_part]
        
        if geometry=='vertex':
            vertices = Part.get_vertices(myPrt, findAt_points, getClosest, searchTolerance)
            myPrt.Set(vertices=vertices, name=name_set)
        
        elif geometry=='edge':
            edges = Part.get_edges(myPrt, findAt_points, getClosest, searchTolerance)
            myPrt.Set(edges=edges, name=name_set)
        
        elif geometry=='face':
            faces = Part.get_faces(myPrt, findAt_points)
            myPrt.Set(faces=faces, name=name_set)
        
        elif geometry=='cell':
            cells = Part.get_cells(myPrt, findAt_points)
            myPrt.Set(cells=cells, name=name_set)
        
    @staticmethod
    def rename_feature(myPrt, new_name, index=-1):
        '''
        Rename a feature object. 
        A feature can be a datum object, partition, reference, etc.
        
        Parameters
        --------------
        new_name: str
            new name of the feature.
            
        index: int
            index of the feature object. 
            Default is -1, represents the last feature created.
        '''
        old_name = myPrt.features.keys()[index]
        
        myPrt.features.changeKey(fromName=old_name, toName=new_name)
        
    @staticmethod
    def create_sketch_point(mySkt, x, y):
        '''
        Create sketch point by coordinates
        
        Parameters
        -------------
        mySkt: 
            Sketch object
            
        x, y: float
            coordinates
        '''
        mySkt.Spot(point=(x, y))
    
    @staticmethod
    def get_datum_by_index(myPrt, index):
        '''
        Get the datum by the index of creation
        
        Parameters
        --------------
        index: int
            the index of datum when it's created (starts from 0).
            
        Returns
        --------------
        datum: Datum object
        '''
        key = myPrt.datums.keys()[index]
        return myPrt.datums[key]

    @staticmethod
    def get_datum_by_name(myPrt, name):
        '''
        Get the datum by name
        
        Parameters
        --------------
        name: str
            name of the datum
            
        Returns
        --------------
        datum: Datum object
        '''
        id_datum = myPrt.features[name].id
        return myPrt.datums[id_datum]

    @staticmethod
    def create_datum_point(myPrt, x, y, z):
        '''
        Create datum point by coordinates
        
        Parameters
        -------------
        myPrt: 
            Part object
            
        x, y, z: float
            coordinates
        '''
        myPrt.DatumPointByCoordinate(coords=(x, y, z))

    @staticmethod
    def create_datum_point_on_edge(myPrt, x, y, z, ratio=0.5):
        '''
        Creates a Feature object and a DatumPoint object along an edge 
        at a selected distance from one end of the edge
        
        Parameters
        -------------
        myPrt: 
            Part object
            
        x, y, z: float
            coordinates to locate the edge
            
        ratio: float
            ratio of the point to one end of the edge
            
        Returns
        -------------
        pt: tuple
            coordinates, (x,y,z)
        '''
        s = myPrt.edges.getClosest(coordinates=((x,y,z),), searchTolerance=100)
        e, closest_point = s[0]
        f = myPrt.DatumPointByEdgeParam(edge=e, parameter=ratio)
        id_datum = f.id
        datum_point = myPrt.datums[id_datum]
        pt = myPrt.getCoordinates(datum_point)
        
        return pt

    @staticmethod
    def create_datum_csys_3p(myPrt, name, origin, dx, dy):
        '''
        Create Datum Coordinate Systems by three points
        
        Parameters
        --------------
        myPrt: 
            Part object
            
        name: str
            name of the datum coordinate system
            
        origin: list, ndarray [3]
            origin coordinates
            
        dx, dy: list, ndarray [3]
            direction vector of the x and y axis
        '''
        origin = np.array(origin)
        dx = np.array(dx)
        dy = np.array(dy)
        
        myPrt.DatumCsysByThreePoints(name=name, coordSysType=CARTESIAN, 
                origin=tuple(origin), point1=tuple(origin + dx), point2=tuple(origin + dy))

    @staticmethod
    def get_vertex(myPrt, findAt_point, toArray=False, getClosest=False, searchTolerance=1000):
        '''
        Get a Vertex/VertexArray object by the `findAt` or `getClosest` command.
        
        -   `findAt` returns the object or objects in the VertexArray located at the given coordinates,
        
        or at a distance of less than 1E-6 from the arbitrary point.
        
        -   `getClosest` returns an object or objects in the VertexArray closest to the given set of points, 
            where the given points need not lie on the Vertex objects in the VertexArray.
        
        http://130.149.89.49:2080/v6.14/books/ker/default.htm?startat=book01.html
        
        Parameters
        -------------
        myPrt: Abaqus Part
            Part object
            
        findAt_point: tuple
            one point coordinate tuple (x,y,z)
            
        toArray: bool
            whether output the Vertex object or VertexArray
            
        getClosest: bool
            whether use `getClosest` to find an vertex
            
        searchTolerance: float
            the distance within which the closest object must lie
            
        Returns
        -------------
        vertex: Abaqus Vertex/VertexArray object
            one vertex

        Note
        -------------
        `myPrt.vertices.findAt((pt))` returns a `Vertex` object.
        
        `myPrt.vertices.findAt((pt,))` returns a `VertexArray Sequence` object.
        
        `myPrt.vertices.findAt((pt,))[0]` is the same as the first one.
        
        `vertices[0]` is a `Vertex` object.
        '''
        if getClosest:
            # The first object in the sequence is an Vertex that is close to the input point referred to by the key. 
            # The second object in the sequence is a sequence of floats that specifies the X-, Y-, and Z-location 
            # of the closest point on the Vertex to the given point. 
            s = myPrt.vertices.getClosest(coordinates=(findAt_point,), searchTolerance=searchTolerance)
            e, closest_point = s[0]
        else:
            # An Vertex object
            e = myPrt.vertices.findAt((findAt_point))
        
        if not toArray:
            return e
        else:
            return myPrt.vertices[e.index:e.index+1]

    @staticmethod
    def get_edge(myPrt, findAt_point, toArray=False, getClosest=False, searchTolerance=1E-6):
        '''
        Get a Edge/EdgeArray object by the `findAt` or `getClosest` command.
        
        -   `findAt` returns the object or objects in the EdgeArray located at the given coordinates,
        
        or at a distance of less than 1E-6 from the arbitrary point.
        
        -   `getClosest` returns an object or objects in the EdgeArray closest to the given set of points, 
            where the given points need not lie on the edges in the EdgeArray.
        
        http://130.149.89.49:2080/v6.14/books/ker/default.htm?startat=pt01ch07pyo03.html
        
        Parameters
        -------------
        myPrt:
            Part object
            
        findAt_point: tuple
            one point coordinate tuple (x,y,z)
            
        toArray: bool
            whether output the Edge object or EdgeArray
            
        getClosest: bool
            whether use `getClosest` to find an edge
            
        searchTolerance: float
            the distance within which the closest object must lie
            
        Returns
        -------------
        edge: Abaqus Edge/EdgeArray object
            one edge

        Note
        -------------
        `myPrt.edges.findAt((pt))` returns an `Edge` object.
        
        `myPrt.edges.findAt((pt,))` returns an `EdgeArray Sequence` object.
        
        `myPrt.edges.findAt((pt,))[0]` is the same as the first one.
        
        `edges[0]` is an `Edge` object.
        '''
        if getClosest:
            # The first object in the sequence is an Edge that is close to the input point referred to by the key. 
            # The second object in the sequence is a sequence of floats that specifies the X-, Y-, and Z-location 
            # of the closest point on the Edge to the given point. 
            s = myPrt.edges.getClosest(coordinates=(findAt_point,), searchTolerance=searchTolerance)
            e, closest_point = s[0]
        else:
            # An Edge object
            e = myPrt.edges.findAt((findAt_point))
        
        if not toArray:
            return e
        else:
            return myPrt.edges[e.index:e.index+1]

    @staticmethod
    def get_face(myPrt, findAt_point, toArray=False):
        '''
        Get a Face/FaceArray object by the `findAt` command.
        
        Parameters
        -------------
        myPrt: Abaqus Part
            Part object
            
        findAt_point: tuple
            one point coordinate tuple (x,y,z)
            
        toArray: bool
            whether output the Face object or FaceArray
            
        Returns
        -------------
        face: Abaqus Face/FaceArray object
            one face

        Note
        -------------
        `myPrt.faces.findAt((pt))` returns a `Face` object.
        
        `myPrt.faces.findAt((pt,))` returns a `FaceArray Sequence` object.
        
        `myPrt.faces.findAt((pt,))[0]` is the same as the first one.
        
        `faces[0]` is a `Face` object.
        '''
        f = myPrt.faces.findAt((findAt_point))
        if not toArray:
            return f
        else:
            return myPrt.faces[f.index:f.index+1]
        
    @staticmethod
    def get_cell(myPrt, findAt_point, toArray=False):
        '''
        Get a Cell/CellArray object by the `findAt` command.
        
        Parameters
        -------------
        myPrt: Abaqus Part
            Part object
        
        findAt_point: tuple
            one point coordinate tuple (x,y,z)
            
        toArray: bool
            whether output the Face object or FaceArray
            
        Returns
        -------------
        cell: Abaqus Cell/CellArray object
            one cell

        Note
        -------------
        `myPrt.cells.findAt((pt))` returns a `Cell Object`.
        
        `myPrt.cells.findAt((pt,))` returns a `CellArray Sequence`.
        
        `myPrt.cells.findAt((pt,))[0]` is the same as the first one.
        
        `cells[0]` is a `Cell` object.
        '''
        c = myPrt.cells.findAt((findAt_point))
        if not toArray:
            return c
        else:
            return myPrt.cells[c.index:c.index+1]

    @staticmethod
    def get_vertices(myPrt, findAt_points, getClosest=False, searchTolerance=1000):
        '''
        Get a VertexArray (Sequence) by the `findAt` or `getClosest` command.
        
        Parameters
        -------------
        myPrt: Abaqus Part
            Part object
        
        findAt_points: tuple, or a list of tuple
        
            point coordinates.
            
            Each vertex is specified by one point coordinate tuple (x,y,z).
            
            `findAt_points` can be either a tuple of one point, or a list of point tuples.
            
        getClosest: bool
            whether use `getClosest` to find an vertex
            
        searchTolerance: float
            the distance within which the closest object must lie
            
        Returns
        -------------
        edges: Abaqus VertexArray (Sequence)
            one or multiple vertex(vertices)
        '''
        if isinstance(findAt_points, tuple):
            findAt_points = [findAt_points]
        
        edges = None
        for pt in findAt_points:
            
            e = Part.get_vertex(myPrt, pt, toArray=True, 
                    getClosest=getClosest, searchTolerance=searchTolerance)
            
            if edges == None:
                edges = e
            else:
                edges += e
        
        return edges

    @staticmethod
    def get_edges(myPrt, findAt_points, getClosest=False, searchTolerance=1000):
        '''
        Get a EdgeArray (Sequence) by the `findAt` or `getClosest` command.
        
        Parameters
        -------------
        myPrt: Abaqus Part
            Part object
            
        findAt_points: tuple, or a list of tuple
        
            point coordinates.
            
            Each edge is specified by one point coordinate tuple (x,y,z).
            
            `findAt_points` can be either a tuple of one point, or a list of point tuples.
            
        getClosest: bool
            whether use `getClosest` to find an edge
            
        searchTolerance: float
            the distance within which the closest object must lie
            
        Returns
        -------------
        edges: Abaqus EdgeArray (Sequence)
            one or multiple edge(s)
        '''
        if isinstance(findAt_points, tuple):
            findAt_points = [findAt_points]
        
        edges = None
        for pt in findAt_points:
            
            e = Part.get_edge(myPrt, pt, toArray=True, 
                    getClosest=getClosest, searchTolerance=searchTolerance)
            
            if edges == None:
                edges = e
            else:
                edges += e
        
        return edges

    @staticmethod
    def get_faces(myPrt, findAt_points):
        '''
        Get a FaceArray (Sequence) by the `findAt` command.
        
        Parameters
        -------------
        myPrt: Abaqus Part
            Part object
            
        findAt_points: tuple, or a list of tuple
        
            point coordinates.
            
            Each face is specified by one point coordinate tuple (x,y,z).
            
            `findAt_points` can be either a tuple of one point, or a list of point tuples.
            
        Returns
        -------------
        faces: Abaqus FaceArray (Sequence)
            one or multiple face(s)
        '''
        if isinstance(findAt_points, tuple):
            findAt_points = [findAt_points]
        
        faces = None
        for pt in findAt_points:
            
            f = Part.get_face(myPrt, pt, toArray=True)
            
            if faces == None:
                faces = f
            else:
                faces += f
        
        return faces

    @staticmethod
    def get_cells(myPrt, findAt_points):
        '''
        Get a CellArray (Sequence) by the `findAt` command.
        
        Parameters
        -------------
        myPrt: Abaqus Part
            Part object
            
        findAt_points: tuple, or a list of tuple
        
            point coordinates.
            
            Each cell is specified by one point coordinate tuple (x,y,z).
            
            `findAt_points` can be either a tuple of one point, or a list of point tuples.
            
        Returns
        -------------
        cells: Abaqus CellArray (Sequence)
            one or multiple cell(s)
        '''
        if isinstance(findAt_points, tuple):
            findAt_points = [findAt_points]
        
        cells = None
        for pt in findAt_points:

            c = Part.get_cell(myPrt, pt, toArray=True)
            
            if cells == None:
                cells = c
            else:
                cells += c
        
        return cells
    
    @staticmethod
    def get_vertex_DatumPointByEdgeParam(myPrt, edge, parameter=0.5):
        '''
        Get the vertex on the edge by `DatumPointByEdgeParam` command.
        
        http://130.149.89.49:2080/v6.14/books/ker/default.htm?startat=pt01ch07pyo04.html
        
        http://130.149.89.49:2080/v6.14/books/ker/default.htm?startat=pt01ch07pyo04.html
        
        Parameters
        --------------
        myPrt: Abaqus Part
            Part object
            
        edge: Abaqus Edge object
            one edge
            
        parameter: float
            specify the distance along edge to the DatumPoint object, in (0,1)
            
        Returns
        ---------------
        point: tuple
        
            coordinate tuple (x,y,z).
            
            Note that this point may not be in myPrt.vertices.
        '''
        myPrt.DatumPointByEdgeParam(edge=edge, parameter=parameter)
        key = myPrt.datums.keys()[-1]
        return myPrt.datums[key].pointOn
    
    @staticmethod
    def get_vertices_on_face(myPrt, face_findAt_point, kind='all', vectors=[]):
        '''
        Get the vertices on the face found at given point.
        
        Parameters
        --------------
        myPrt: Abaqus Part
            Part object
            
        findAt_points: tuple
            point coordinate tuple (x,y,z), the face is specified by one point.
            
        kind: str
        
            specify vertices to be obtained.

            - 'all': all vertices.
            - 'vec': use vectors to find vertices.
            - 'x*y*z*': vertex specified by minimum (I) or maximum (A) x,y,z coordinates. '*' is 'I' or 'A' or 'N'.
        
        vectors: list
            vectors to locate each point by finding the furthest point in that direction, [[3],[3],...]
            
        Returns
        ---------------
        vertices: list of tuple
            a list of vertex coordinate tuple (x,y,z)
        
        Attributes
        ---------------
        id_vertices: list
            vertex ids of the vertices of the face
            
        all_vertices: list
            a list of coordinate tuple (x,y,z) for all vertices
        '''
        #* Get Face Object
        face = Part.get_face(myPrt, face_findAt_point)

        #* All vertices on the face
        id_vertices = face.getVertices()
        all_vertices = []
        for id in id_vertices:
            all_vertices.append(myPrt.vertices[id].pointOn[0])
    
        #* Pick vertices
        if kind == 'all':
            return all_vertices
        
        if kind == 'vec':
            
            if len(vectors)==0:
                raise Exception
            
            vertices = []
            for project_vec in vectors:
                proj_len = np.dot(vertices, project_vec)
                ii = np.argmax(proj_len)
                vertices.append(all_vertices[ii])
                
            return vertices
        
        if not isinstance(kind, list):
            kind = [kind]
        
        vertices = []
        for kd in kind:
            v, _ = Part.get_vertex_xSySzS(all_vertices, kind=kd)
            vertices.append(v)
            
        return vertices
    
    @staticmethod 
    def get_vertex_xSySzS(vertices, kind='xIyIzI'):
        '''
        Get the vertex in a list of vertices by minimum/maximum x,y,z.
        
        Parameters
        --------------
        vertices: list
            a list of coordinate tuple (x,y,z) for all vertices
            
        kind: str
        
            specify vertices to be obtained.
            
            'x*y*z*': vertex specified by minimum (I) or maximum (A) x,y,z coordinates. 
            '*' is 'I' or 'A' or 'N'.

        Returns
        ---------------
        vertex: tuple
            coordinate tuple (x,y,z)
            
        index: int
            index of the vertex in vertices
        '''
        vertices = np.array(vertices)
        project_vec = np.zeros(3)
        
        for i in range(3):
            if kind[2*i+1] == 'A':
                project_vec[i] = 1.0
            elif kind[2*i+1] == 'I':
                project_vec[i] = -1.0
            else:
                project_vec[i] = 0.0
                
        proj_len = np.dot(vertices, project_vec)
        ii = np.argmax(proj_len)
        
        return vertices[ii], ii
    
    #* =============================================
    #* Abaqus Mesh functions
    def get_CompositeLayup_thickness(self, name_set, total_thickness, ply_angle, eNum_thickness=1, symmetric=True):
        '''
        Get the list of thickness of each ply in the composite layup.
        
        Parameters
        ------------
        name_set: str
            name of the set (cells) for Error information
            
        total_thickness: float
            total thickness of the composite plate
            
        ply_angle: list
            orientation angle of each ply
        
        eNum_thickness: int
            number of elements in the thickness direction
        
        symmetric: bool
            whether the composite is symmetric
        
        Returns
        ------------
        thickness: list
            thickness of each ply
        '''
        ply_thickness = self.pMesh['composite_ply_thickness']
        num_ply = len(ply_angle)
        ratio = 2 if symmetric else 1
        
        if ply_thickness <= 0:
            
            t = total_thickness / float(num_ply*ratio*eNum_thickness)
            thickness = [t for _ in range(num_ply)]
            
            print('Automatically set ply thickness = %.3f'%(t))
            
        else:
        
            thickness = [ply_thickness for _ in range(num_ply)]

            valid = np.sum(thickness)*ratio*eNum_thickness == total_thickness
        
            if not valid and self.raise_exception:
                print('Error: Set [%s]'%(name_set))
                print('       The composite plate thickness does not match the number and thickness of plies')
                print('       symmetric laminate            = %d'%(int(symmetric)))
                print('       thickness direction elements  = %d'%(int(symmetric)))
                print('       number of plies               = %d'%(int(num_ply*ratio)))
                print('       total thickness of laminate   = %.3f'%(total_thickness))
                print('       sum(ply thickness)            = %.3f'%(np.sum(thickness)*ratio*eNum_thickness))
                raise Exception()
            
        if np.min(thickness) <= 0 and self.raise_exception:
            print('Error: Set [%s]'%(name_set))
            print('       The ply thickness should be positive')
            print(thickness, total_thickness, ply_thickness, num_ply, ratio, eNum_thickness)
            raise Exception()
        
        return thickness

    def set_CompositeLayup_of_set(self, myPrt, name_set, total_thickness, ply_angle, 
                    eNum_thickness=1, symmetric=True, numIntPoints=3, name_csys_datum=None, 
                    stackDirection=None, normalDirection=None,
                    rotation_angle=None, material='IM7/8551-7'):
        '''
        Edit Composite Layup for a set.
        
        Parameters
        --------------------
        myPrt:
            Abaqus Part object
            
        name_set: str
            name of the set
            
        total_thickness: float
            thickness of the plate
            
        ply_angle: list
            angles of each ply (degree)
            
        eNum_thickness: int
            number of elements in the thickness direction
            
        symmetric: bool
            whether the layup is symmetric
            
        numIntPoints: int
            number of integration points
            
        name_csys_datum: str or None
            name of the datum coordinate system for the layup orientation
            
        stackDirection: Abaqus constant or None
        
            `STACK_3`: default (when the input is None), 
            it stands for 'element direction 3', which is the bottom-to-top direction 
            defined in `Mesh -> Assign Stack Direction`.
            
            `STACK_ORIENTATION`: it stand for 'layup orientation'
            
        normalDirection: Abaqus constant or None

            `AXIS_3`: default (when the input is None),
            it stands for the 3rd axis (z axis) of the datum or global coordinate system.
            
            `AXIS_1`, `AXIS_2`: the 1st, 2nd axis of the datum or global coordinate system.
            
        rotation_angle: float or None
            the rotation angle of the datum coordinate system about the normal axis (3) 
            for the layup orientation
            
        material: str
            name of the composite material. It is 'IM7/8551-7' by default.
        '''
        if name_set not in myPrt.sets.keys():
            return        
        
        # CompositeLayup object
        # https://docs.software.vt.edu/abaqusv2022/English/?show=SIMACAEKERRefMap/simaker-c-compositelayuppyc.htm
        
        compositeLayup = myPrt.CompositeLayup(name=name_set, 
                        description='', elementType=CONTINUUM_SHELL, symmetric=symmetric)
        
        # CompositeShellSection object
        # https://docs.software.vt.edu/abaqusv2022/English/?show=SIMACAEKERRefMap/simaker-c-compositeshellsectioncpp.htm
        #
        # The default poisson ratio is 0.5
        # The thicknessModulus specifies the effective thickness modulus, it is relevant only for continuum shells.
        # The default value is twice the initial in-plane shear modulus based on the material definition.
        # 
        # # https://docs.software.vt.edu/abaqusv2022/English/SIMACAEELMRefMap/simaelm-c-shellelem.htm#simaelm-c-shellelem-contshellthick
        #
        # The default thickness modulus is computed as the thickness-weighted harmonic mean of the contributions from individual layers.
        # The thickness modulus of a layer is twice the initial in-plane elastic shear modulus based on the material definition 
        # for that layer in the initial configuration.
        #
        compositeLayup.Section(preIntegrate=OFF, 
                        integrationRule=SIMPSON, poissonDefinition=DEFAULT, 
                        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF)
    
        # MaterialOrientation object
        # https://docs.software.vt.edu/abaqusv2022/English/?show=SIMACAEKERRefMap/simaker-c-materialorientationpyc.htm
    
        #* Layup Orientation Attributes
        if stackDirection is None:
            stackDirection = STACK_3
            
        if normalDirection is None:
            normalDirection = AXIS_3
                
        #* Layup Orientation: datum coordinate system
        if name_csys_datum is not None:
            
            csys_datum = self.get_datum_by_name(myPrt, name_csys_datum)

            compositeLayup.orientation.setValues(
                orientationType=SYSTEM,             # Definition: coordinate system
                localCsys=csys_datum,               # Select CSYS
                axis=normalDirection,               # Normal direction: axis 3
                stackDirection=stackDirection)      # Stacking direction: layup orientation
            
            if rotation_angle is not None:
                compositeLayup.orientation.setValues(
                    additionalRotationType=ROTATION_ANGLE, 
                    angle=rotation_angle, additionalRotationField='')
        
        #* Layup Orientation: part global
        else:

            compositeLayup.ReferenceOrientation(
                orientationType=GLOBAL, localCsys=None, fieldName='', 
                additionalRotationType=ROTATION_NONE, angle=0.0, 
                axis=normalDirection,           # Normal direction: part global axis **
                stackDirection=STACK_3)         # Stacking direction: element direction 3

        #* Plies
        thickness = self.get_CompositeLayup_thickness(
                name_set, total_thickness, ply_angle, eNum_thickness, symmetric)
        
        for i in range(len(ply_angle)):
        
            compositeLayup.CompositePly(suppressed=False, plyName='ply-%d'%(i+1), 
                region=myPrt.sets[name_set], material=material, 
                thicknessType=SPECIFY_THICKNESS, thickness=thickness[i], 
                orientationType=SPECIFY_ORIENT, orientationValue=ply_angle[i], 
                additionalRotationType=ROTATION_NONE,   # CSYS in `Plies` for rotation angle
                additionalRotationField='',             # CSYS in `Plies` for rotation angle
                axis=AXIS_3, angle=0.0,                 # CSYS in `Plies` for rotation angle
                numIntPoints=numIntPoints)

    @staticmethod
    def get_seedEdgeByBias_ratio(min_size, edge_length, n_element, biasMethod='single', factor=1.5):
        '''
        Calculate the `ratio` for `seedEdgeByBias(ratio=ratio, ***)`.
        
        Parameters
        --------------
        min_size: float
            minimum element size
            
        edge_length: float
            total length of the edge
            
        n_element: int
            number of element
            
        biasMethod: str
            'single' or 'double'
            
        factor: float
            multiply the ratio by `factor`
        
        Returns
        --------------
        ratio: float
            A Float specifying the ratio of the largest element to the smallest element. 
            Possible values are 1.0 ~ 1E6.
        '''
        ratio = 1.0
        
        if biasMethod == 'single':
            aa = edge_length/(min_size*1.0)
        else:
            aa = edge_length/(min_size*2.0)
                
        def func(x):
            return aa*(x-1) - (x**n_element-1)
        
        for q0 in [1.2, 1.5, 2.0, 5.0, 10.0]:
            if func(q0)<=0.0:
                break
        
        q = fsolve(func, q0)[0]
        ratio = q**(n_element-1)
        ratio = np.ceil(10*factor*ratio)*0.1
        
        if q<1.0:
            print('[Error] Calculating ratio for seedEdgeByBias(ratio=ratio, ***)')
            print('        q should be larger than 1, yet q = ', q)
            print('        ratio = ', ratio)
            raise Exception
        
        return ratio

    @staticmethod
    def set_element_type_of_part(myPrt, name_set=None, kind='continuum shell', elemLibrary='standard', 
                        hourglassControl='default', maxDegradation=None):
        '''
        Set element type for a part.
        
        Parameters
        --------------
        name_set: str | None
            name of the set for the element type, or None for the whole part
            
        kind: str
            'continuum shell', '3D stress', 'Cohesive'
            
        elemLibrary: str
            'standard' for static simulation, 
            'explicit' for dynamic simulation
        
        hourglassControl:  str
            'default' for 'DEFAULT',
            'enhanced' for 'ENHANCED'
        
        maxDegradation: float or None
            max degradation (Dmax) controls how Abaqus treats elements with severe damage. 
            Dmax is the upper bound to the overall damage variable D; 
            and you can choose whether to delete an element once maximum degradation is reached. 
            The latter choice also affects which stiffness components are damaged.
            
            For example, Dmax = 0.99.
            
            The default setting of Dmax depends on whether elements are to be deleted upon reaching maximum degradation (discussed next). 
            For the default case of element deletion and in all cases for cohesive elements, Dmax=1.0; otherwise, Dmax=0.99. 
            The output variable SDEG contains the value of D. No further damage is accumulated at an integration point once D reaches Dmax 
            (except, of course, any remaining stiffness is lost upon element deletion).
        '''
        if elemLibrary=='standard':
            elemLibrary = STANDARD
        elif elemLibrary=='explicit':
            elemLibrary = EXPLICIT
        else:
            raise Exception
        
        if hourglassControl=='default':
            hourglassControl = DEFAULT
        elif hourglassControl=='enhanced':
            hourglassControl = ENHANCED
        else:
            raise Exception
        
        if kind=='continuum shell':
            
            if maxDegradation is None:
                elemType1 = mesh.ElemType(elemCode=SC8R, elemLibrary=elemLibrary, 
                                            secondOrderAccuracy=OFF, hourglassControl=hourglassControl)
            else:
                elemType1 = mesh.ElemType(elemCode=SC8R, elemLibrary=elemLibrary, 
                                            secondOrderAccuracy=OFF, hourglassControl=hourglassControl, 
                                            maxDegradation=maxDegradation)
                
            elemType2 = mesh.ElemType(elemCode=SC6R, elemLibrary=elemLibrary)
            elemType3 = mesh.ElemType(elemCode=UNKNOWN_TET, elemLibrary=elemLibrary)
            
        elif kind=='3D stress': # (solid cells)
            
            elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=elemLibrary, 
                                        kinematicSplit=AVERAGE_STRAIN, 
                                        secondOrderAccuracy=OFF, hourglassControl=hourglassControl, 
                                        distortionControl=DEFAULT)
            elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=elemLibrary)
            elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=elemLibrary)
            
        elif kind=='Cohesive': # (cohesive cells)
            
            elemType1 = mesh.ElemType(elemCode=COH3D8, elemLibrary=STANDARD, viscosity=1e-05)
            elemType2 = mesh.ElemType(elemCode=COH3D6, elemLibrary=STANDARD)
            elemType3 = mesh.ElemType(elemCode=UNKNOWN_TET, elemLibrary=STANDARD)
            
        else:
            
            raise Exception

        # Assign element types
        if name_set is None:
            regions = (myPrt.cells,)
        else:
            regions = (myPrt.sets[name_set].cells,)
        
        myPrt.setElementType(regions=regions, 
                                elemTypes=(elemType1, elemType2, elemType3))

