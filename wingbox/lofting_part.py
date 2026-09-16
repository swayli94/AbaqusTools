'''

'''
import os
import time
import numpy as np

from AbaqusTools import Part, IS_ABAQUS

if IS_ABAQUS:
    from abaqus import *
    from abaqusConstants import *
    from sketch import *
    import mesh

path = os.path.dirname(os.path.abspath(__file__))

from geometry import WingSectionGeometry, get_primaryAxisVector_spanwise
from utils import mid_pt, mid_pts
from layup import create_shell_CompositeLayup_of_set


def get_component_ply_thickness(pMesh, params):
    '''
    Get component-level ply thickness, with old global value as fallback.
    '''
    return params.get('ply_thickness', pMesh['ply_thickness'])


def is_non_design_wingbox_section(pMesh, i_section):
    '''
    Return True when a spanwise wingbox bay is part of the non-design region.
    '''
    return int(i_section) in [int(i) for i in pMesh.get('non_design_wingbox_sections', [])]


def get_non_design_thickness_factor(pMesh, component, i_section):
    '''
    Get the ply-thickness multiplier for a non-design wingbox bay.
    '''
    if not is_non_design_wingbox_section(pMesh, i_section):
        return 1.0

    non_design_region = pMesh.get('non_design_region', {})
    factors = non_design_region.get('thickness_factor', {})
    if component in factors:
        return float(factors[component])
    if 'default' in factors:
        return float(factors['default'])
    return 1.0


def get_design_region_ply_thickness(pMesh, params, component, i_section):
    '''
    Get component ply thickness with optional non-design region thickening.
    '''
    thickness = get_component_ply_thickness(pMesh, params)
    return thickness * get_non_design_thickness_factor(pMesh, component, i_section)


def get_cover_layup_params(pMesh, side):
    '''
    Get upper/lower cover layup parameters.

    Supports both:
    - old format: pMesh['cover']['layup_orientAngles']
    - new format: pMesh['cover']['upper'/'lower']['layup_orientAngles']
    '''
    cover = pMesh['cover']
    if side in cover:
        return cover[side]
    return cover


def get_span_group_layup_params(pMesh, params, i_section):
    '''
    Resolve a component's layup parameters to the span group of a bay.

    Supports both:
    - old format: no 'span_groups' in pMesh; `params` is the layup and is
      returned unchanged, so parameter files written before the spanwise
      split keep working;
    - span-group format: `pMesh['span_groups']` maps each group name to
      the bays it covers and `params` maps the same names to the group's
      layup; the layup of the group covering bay `i_section` is returned.

    A bay no group covers raises ValueError: it would otherwise be built
    from a laminate that was never chosen for it.
    '''
    span_groups = pMesh.get('span_groups')
    if span_groups is None:
        return params
    for group, bays in span_groups.items():
        if int(i_section) in bays:
            return params[group]
    raise ValueError(
        'No span group covers bay %d (span_groups: %s).'
        % (int(i_section), span_groups)
    )


class LoftingPart(Part):
    '''
    Class for wingbox lofting part.
    '''
    def __init__(self, name_part, model, pGeo, pMesh):
        super(LoftingPart,self).__init__(model, pGeo, pMesh)
        self.name_part = name_part

        self.sections = []
        for section_params in pGeo['sections']:
            section_params['airfoil'] = os.path.join(path, section_params['airfoil'])
            wsg = WingSectionGeometry()
            wsg.set_parameters(section_params)
            self.sections.append(wsg)
            
        self.name_layups = []

    def create_sketch(self):
        '''
        Create stand-alone sketches via Abaqus Module: Sketch
        
        - Create Sketch
        '''
        self._create_sketch_section_for_sweep()
        self._create_sketch_reference_plane()

    def _create_sketch_section_for_sweep(self):
        '''
        Create stand-alone sketch for sweeping: cover + spar + stringer
        '''
        #* Create the open-hole plate sketch in x-y plane
        for i, section in enumerate(self.sections):
            
            mySkt = self.model.ConstrainedSketch(name='section_%d'%(i), sheetSize=section.chord*3)
            
            cover_u = np.concatenate((section.x3d_upper_cover[:,np.newaxis], section.y3d_upper_cover[:,np.newaxis]), axis=1) # [n,2]
            cover_l = np.concatenate((section.x3d_lower_cover[:,np.newaxis], section.y3d_lower_cover[:,np.newaxis]), axis=1) # [n,2]
            
            mySkt.Spline(points=tuple(cover_u))
            mySkt.Spline(points=tuple(cover_l))
            
            for j, spar in enumerate(section.spars):
                pt1 = spar.get_selection_point(feature='root', side='upper')[:2]
                pt2 = spar.get_selection_point(feature='root', side='lower')[:2]
                mySkt.Line(point1=pt1, point2=pt2)
                
            for j, stringer in enumerate(section.stringers):
                for side in ['upper', 'lower']:
                    pt1 = stringer.get_selection_point(feature='root', side=side)[:2]
                    pt2 = stringer.get_selection_point(feature='corner', side=side)[:2]
                    pt3 = stringer.get_selection_point(feature='tip', side=side)[:2]
                    mySkt.Line(point1=pt1, point2=pt2)
                    mySkt.Line(point1=pt2, point2=pt3)

    def _create_sketch_reference_plane(self):
        '''
        Create a stand-alone sketch for reference plane,
        used for shell lofting.
        '''
        mySkt = self.model.ConstrainedSketch(name='ref', sheetSize=1000.0)
        mySkt.rectangle(point1=(0.0, 0.0), point2=(10.0, 10.0))
        
    def create_part(self):
        '''
        Create part by shell lofting
        '''
        myPrt = self.model.Part(name=self.name_part, dimensionality=THREE_D, type=DEFORMABLE_BODY)

        #* Reference plane and axis
        myPrt.DatumAxisByPrincipalAxis(principalAxis=XAXIS)
        self.rename_feature(myPrt, 'XAXIS')
        myPrt.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
        self.rename_feature(myPrt, 'YAXIS')
        myPrt.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
        self.rename_feature(myPrt, 'ZAXIS')
        self.create_datum_csys_3p(myPrt, 'csys_plate', origin=[0.0, 0.0, 0.0],
                                    dx=[1, 0, 0], dy=[0, 1, 0])
        
        for i, section in enumerate(self.sections):
            myPrt.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=section.zLE)
            self.rename_feature(myPrt, 'XYPLANE-%d'%(i))

        #* Plane for plate sketch
        transform = myPrt.MakeSketchTransform(
            sketchPlane=self.get_datum_by_name(myPrt, 'XYPLANE-0'),
            sketchUpEdge=self.get_datum_by_name(myPrt, 'XAXIS'), 
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 0.0))
    
        #* Section sketch
        max_chord = max([section.chord for section in self.sections])
        mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=max_chord*10, transform=transform)
        mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)
        mySkt.retrieveSketch(sketch=self.model.sketches['ref'])
        
        #* Create part by Shell (Planer)
        myPrt.BaseShell(sketch=mySkt)
        self.rename_feature(myPrt, 'reference_plane')
        del self.model.sketches['__profile__']
        
        #* ============================================
        #* Create lofting for each two sections
        #* ============================================
        n_sections = len(self.sections)
        
        #* Add lofting sections as wires
        for i in range(n_sections):
            section = self.sections[i]

            transform = myPrt.MakeSketchTransform(
                sketchPlane=self.get_datum_by_name(myPrt, 'XYPLANE-%d'%(i)),
                sketchUpEdge=self.get_datum_by_name(myPrt, 'XAXIS'), 
                sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 0.0))
        
            mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=section.chord*3, transform=transform)
            mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)
            mySkt.retrieveSketch(sketch=self.model.sketches['section_%d'%(i)])
            
            myPrt.Wire(sketchPlane=self.get_datum_by_name(myPrt, 'XYPLANE-%d'%(i)), 
                sketchUpEdge=self.get_datum_by_name(myPrt, 'XAXIS'), 
                sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, sketch=mySkt)
            self.rename_feature(myPrt, 'wire_section_%d'%(i))
            del self.model.sketches['__profile__']
        
        #* Lofting for each two sections
        for i_section in range(n_sections-1):
            
            sections = [self.sections[i_section], self.sections[i_section+1]]

            #* Lofting for covers
            for side in ['upper', 'lower']:
                loftsections = [] # [(edge1, edge2, ...), (edge3, edge4, ...) ...]
                for i, section in enumerate(sections):
                    findAt_points = section.get_selection_points(feature='cover', side=side, index=0)
                    try:
                        # Arbitrary stringer x/c positions are interpolated on
                        # the analytical airfoil coordinates, while Abaqus
                        # creates a spline through the airfoil samples.  Their
                        # points can differ slightly, so strict findAt (1e-6
                        # model units) is not reliable for refined layouts.
                        # One millimetre is still tiny compared with the
                        # minimum reinforcement-line spacing (>300 mm here).
                        edges = self.get_edges(
                            myPrt,
                            findAt_points=findAt_points,
                            getClosest=True,
                            searchTolerance=1.0)
                    except Exception as exc:
                        raise RuntimeError(
                            'Cover wire-edge lookup failed: '
                            'section_pair=%d endpoint=%d side=%s points=%s; '
                            'cause=%s' %
                            (i_section, i, side, findAt_points, exc))
                    loftsections.append(tuple(edges))
                myPrt.ShellLoft(loftsections=tuple(loftsections), startCondition=NONE, endCondition=NONE)
                self.rename_feature(myPrt, 'wingbox%d_cover_%s'%(i_section, side))
            
            #* Lofting for spars
            n_spars = sections[0].n_spars
            for j in range(n_spars):
                loftsections = [] # [(edge1, edge2, ...), (edge3, edge4, ...) ...]
                for i, section in enumerate(sections):
                    findAt_points = section.get_selection_points(feature='spar', side=None, index=j)
                    edges = self.get_edges(myPrt, findAt_points=findAt_points)
                    loftsections.append(tuple(edges))
                myPrt.ShellLoft(loftsections=tuple(loftsections), startCondition=NONE, endCondition=NONE)
                self.rename_feature(myPrt, 'wingbox%d_spar_%d'%(i_section, j))
            
            #* Lofting for stringers
            n_stringers = sections[0].n_stringers
            for j in range(n_stringers):
                for side in ['upper', 'lower']:
                    loftsections = [] # [(edge1, edge2, ...), (edge3, edge4, ...) ...]
                    for i, section in enumerate(sections):
                        findAt_points = section.get_selection_points(feature='stringer', side=side, index=j)
                        try:
                            edges = self.get_edges(
                                myPrt, findAt_points=findAt_points)
                        except Exception as exc:
                            raise RuntimeError(
                                'Stringer wire-edge lookup failed: '
                                'section_pair=%d endpoint=%d stringer_index=%d '
                                'side=%s x_relative=%s points=%s; cause=%s'
                                % (i_section, i, j, side,
                                   getattr(section.stringers[j], 'x', None),
                                   findAt_points, exc))
                        loftsections.append(tuple(edges))
                    myPrt.ShellLoft(loftsections=tuple(loftsections), startCondition=NONE, endCondition=NONE)
                    self.rename_feature(myPrt, 'wingbox%d_stringer_%d_%s'%(i_section, j, side))
        
        #* Adjacent web faces of one spar or stringer are exactly coplanar
        #* when the planform has no twist gradient and a constant sweep,
        #* and the loft kernel then keeps no junction edge at the shared
        #* station.  The seeding and tie lookups in `create_set` need that
        #* edge.  A web face that answers the mid-bay point of both
        #* adjacent bays spans the station; partition it with the
        #* station's datum plane.  Faces the kernel already split are
        #* left alone, so the topology does not depend on the twist
        #* gradient happening to be nonzero.
        for i in range(1, n_sections - 1):
            faces = []
            for j in range(self.sections[i].n_spars):
                faces += self._faces_spanning_station(myPrt, 'spar', j, i)
            for j in range(self.sections[i].n_stringers):
                for side in ['upper', 'lower']:
                    faces += self._faces_spanning_station(
                        myPrt, 'stringer', j, i, side=side)
            if faces:
                myPrt.PartitionFaceByDatumPlane(
                    datumPlane=self.get_datum_by_name(myPrt, 'XYPLANE-%d' % i),
                    faces=faces)
        
        #* Delete the reference plane
        myPrt.setValues(geometryRefinement=EXTRA_FINE)
        del myPrt.features['reference_plane']

    def _faces_spanning_station(self, myPrt, feature, index, i, side=None):
        '''
        Web faces of feature `index` that span the station plane at
        section `i`, i.e. that answer the mid-bay lookup of both adjacent
        bays.  The face recipes are the ones `create_set` uses for the
        per-bay face sets.
        '''
        face_search_tolerance = self.pMesh.get('face_search_tolerance', 1E-2)
        pts_left = mid_pts(
            self.sections[i-1].get_selection_points(feature=feature, side=side, index=index),
            self.sections[i].get_selection_points(feature=feature, side=side, index=index))
        pts_right = mid_pts(
            self.sections[i].get_selection_points(feature=feature, side=side, index=index),
            self.sections[i+1].get_selection_points(feature=feature, side=side, index=index))
        faces_left = self.get_faces(myPrt, pts_left,
            getClosest=True, searchTolerance=face_search_tolerance)
        faces_right = self.get_faces(myPrt, pts_right,
            getClosest=True, searchTolerance=face_search_tolerance)
        right_indices = set(f.index for f in faces_right)
        return [f for f in faces_left if f.index in right_indices]

    def create_surface(self):
        '''
        Create surfaces, including:
        - cover surfaces (upper/lower)
        - spar surfaces (each spar)
        - stringer surfaces (each stringer, upper/lower)
        '''
        n_sections = len(self.sections)
        face_search_tolerance = self.pMesh.get('face_search_tolerance', 1E-2)

        myPrt = self.model.parts[self.name_part]

        for i_section in range(n_sections - 1):
            sec0 = self.sections[i_section]
            sec1 = self.sections[i_section + 1]
            tag = 'wingbox%d' % i_section

            # Cover faces (one face per cover segment between spars/stringers)
            for side in ['upper', 'lower']:
                pts0 = sec0.get_selection_points(feature='cover', side=side, index=0)
                pts1 = sec1.get_selection_points(feature='cover', side=side, index=0)
                findAt_points=mid_pts(pts0, pts1)
                faces = self.get_faces(myPrt, findAt_points,
                    getClosest=True, searchTolerance=face_search_tolerance)
                myPrt.Surface(side1Faces=faces, name='face_%s_cover_%s' % (tag, side))

            # Spar faces
            for j in range(sec0.n_spars):
                pts0 = sec0.get_selection_points(feature='spar', side=None, index=j)
                pts1 = sec1.get_selection_points(feature='spar', side=None, index=j)
                findAt_points=mid_pts(pts0, pts1)
                faces = self.get_faces(myPrt, findAt_points,
                    getClosest=True, searchTolerance=face_search_tolerance)
                myPrt.Surface(side1Faces=faces, name='face_%s_spar_%d' % (tag, j))

            # Stringer faces (web + flange per side, combined into one set)
            for j in range(sec0.n_stringers):
                for side in ['upper', 'lower']:
                    pts0 = sec0.get_selection_points(feature='stringer', side=side, index=j)
                    pts1 = sec1.get_selection_points(feature='stringer', side=side, index=j)
                    findAt_points=mid_pts(pts0, pts1)
                    faces = self.get_faces(myPrt, findAt_points,
                        getClosest=True, searchTolerance=face_search_tolerance)
                    myPrt.Surface(side1Faces=faces, name='face_%s_stringer_%d_%s' % (tag, j, side))

    def create_set(self):
        '''
        Create sets for all segments of the lofting_part, including:
        - cover faces (upper/lower)
        - spar faces (each spar)
        - stringer web, flange faces (each stringer, upper/lower)
        - spar edges (upper/lower edge of each spar)
        - stringer edges (root/corner/tip edge of each stringer)
        '''
        n_sections = len(self.sections)
        face_search_tolerance = self.pMesh.get('face_search_tolerance', 1E-2)
        cover_edge_search_tolerance = self.pMesh.get(
            'cover_edge_search_tolerance', 1.0)

        # ============================================================
        # Faces and edges for each wingbox segment
        # For seeding and section assignment
        # ============================================================
        myPrt = self.model.parts[self.name_part]
        myPrt.Set(faces=myPrt.faces, name='all') 

        for i_section in range(n_sections - 1):
            sec0 = self.sections[i_section]
            sec1 = self.sections[i_section + 1]
            tag = 'wingbox%d' % i_section

            # Cover faces (one face per cover segment between spars/stringers)
            for side in ['upper', 'lower']:
                pts0 = sec0.get_selection_points(feature='cover', side=side, index=0)
                pts1 = sec1.get_selection_points(feature='cover', side=side, index=0)
                findAt_points=mid_pts(pts0, pts1)

                self.create_geometry_set(
                    name_set='face_%s_cover_%s' % (tag, side),
                    findAt_points=findAt_points,
                    geometry='face',
                    getClosest=True, searchTolerance=face_search_tolerance)

            # Spar faces
            for j in range(sec0.n_spars):
                pts0 = sec0.get_selection_points(feature='spar', side=None, index=j)
                pts1 = sec1.get_selection_points(feature='spar', side=None, index=j)
                self.create_geometry_set(
                    name_set='face_%s_spar_%d' % (tag, j),
                    findAt_points=mid_pts(pts0, pts1),
                    geometry='face',
                    getClosest=True, searchTolerance=face_search_tolerance)

            # Stringer faces (web + flange per side, combined into one set)
            for j in range(sec0.n_stringers):
                for side in ['upper', 'lower']:
                    pts0 = sec0.get_selection_points(feature='stringer', side=side, index=j)
                    pts1 = sec1.get_selection_points(feature='stringer', side=side, index=j)
                    self.create_geometry_set(
                        name_set='face_%s_stringer_%d_%s' % (tag, j, side),
                        findAt_points=mid_pts(pts0, pts1),
                        geometry='face',
                        getClosest=True, searchTolerance=face_search_tolerance)

            # Spar edges (upper/lower junction with cover, longitudinal in z)
            for j in range(sec0.n_spars):
                for side in ['upper', 'lower']:
                    pt0 = sec0.spars[j].get_selection_point(feature='root', side=side)
                    pt1 = sec1.spars[j].get_selection_point(feature='root', side=side)
                    self.create_geometry_set(
                        name_set='edge_%s_spar_%d_%s' % (tag, j, side),
                        findAt_points=[mid_pt(pt0, pt1)],
                        geometry='edge')

            # Stringer edges (root/corner/tip, upper/lower, longitudinal in z)
            for j in range(sec0.n_stringers):
                for side in ['upper', 'lower']:
                    for feat in ['root', 'corner', 'tip']:
                        pt0 = sec0.stringers[j].get_selection_point(feature=feat, side=side)
                        pt1 = sec1.stringers[j].get_selection_point(feature=feat, side=side)
                        self.create_geometry_set(
                            name_set='edge_%s_stringer_%d_%s_%s' % (tag, j, side, feat),
                            findAt_points=[mid_pt(pt0, pt1)],
                            geometry='edge')

        # ============================================================
        # Side edges (in xy-plane) for each wingbox section
        # For seeding and tie constraints
        # ============================================================
        for i_section, sec in enumerate(self.sections):
            tag = 'sec%d' % i_section
            
            # Cover side edges
            for side in ['upper', 'lower']:
                pts0 = sec.get_selection_points(feature='cover', side=side, index=0)
                self.create_geometry_set(
                    name_set='edge_%s_cover_%s' % (tag, side),
                    findAt_points=pts0,
                    geometry='edge', getClosest=True,
                    searchTolerance=cover_edge_search_tolerance)
            
            # Spar side edges
            for j in range(sec.n_spars):
                for side in ['upper', 'lower']:
                    pt0 = sec.spars[j].get_selection_point(feature='spar', side=None)
                    self.create_geometry_set(
                        name_set='edge_%s_spar_%d' % (tag, j),
                        findAt_points=pt0,
                        geometry='edge')
            
            # Stringer side edges
            for j in range(sec.n_stringers):
                for side in ['upper', 'lower']:
                    for feat in ['web', 'flange']:
                        pt0 = sec.stringers[j].get_selection_point(feature=feat, side=side)
                        self.create_geometry_set(
                            name_set='edge_%s_stringer_%d_%s_%s' % (tag, j, side, feat),
                            findAt_points=pt0,
                            geometry='edge')

    def set_seeding(self):
        
        myPrt = self.model.parts[self.name_part]
        myPrt.seedPart(size=self.pMesh['global_seedPart_size'], 
                        deviationFactor=0.1, minSizeFactor=0.1)
        
        # ============================================================
        # side edges (in xy-plane) for each wingbox section
        # ============================================================
        for i_section, sec in enumerate(self.sections):
            tag = 'sec%d' % i_section
            
            # Spar side edges (for seeding)
            for j in range(sec.n_spars):
                for side in ['upper', 'lower']:
                    name_set = 'edge_%s_spar_%d' % (tag, j)
                    myPrt.seedEdgeByNumber(edges=myPrt.sets[name_set].edges,
                                    number=self.pMesh['spar_seedEdge_number'], constraint=FIXED)
                    
            # Stringer side edges (for seeding)
            for j in range(sec.n_stringers):
                for side in ['upper', 'lower']:
                    for feat in ['web', 'flange']:
                        name_set = 'edge_%s_stringer_%d_%s_%s' % (tag, j, side, feat)
                        myPrt.seedEdgeByNumber(edges=myPrt.sets[name_set].edges,
                                    number=self.pMesh['stringer_seedEdge_number'], constraint=FIXED)

        self._ignore_redundant_short_edge_vertices(myPrt)

    def _ignore_redundant_short_edge_vertices(self, myPrt):
        """Merge only short, collinear edge fragments using virtual topology.

        A candidate vertex must join exactly two edges, the two edges must have
        the same pair of adjacent faces, and their tangents at the vertex must
        be nearly opposite.  Junction vertices (for example cover/stringer
        attachment points) are therefore never ignored.
        """
        if not bool(self.pMesh.get('short_edge_virtual_topology', False)):
            return

        threshold = float(self.pMesh.get('short_edge_threshold', 1.0))
        cosine_limit = float(
            self.pMesh.get('short_edge_collinearity_cosine', -0.99))
        if threshold <= 0.0:
            raise ValueError('short_edge_threshold must be positive.')

        candidates = {}
        short_edge_records = []
        for short_edge in myPrt.edges:
            short_size = float(short_edge.getSize(printResults=False))
            if short_size >= threshold:
                continue
            short_edge_records.append((int(short_edge.index), short_size))
            short_faces = set([int(i) for i in short_edge.getFaces()])
            if len(short_faces) != 2:
                continue

            for vertex_index in short_edge.getVertices():
                vertex_index = int(vertex_index)
                vertex = myPrt.vertices[vertex_index]
                incident_indices = [int(i) for i in vertex.getEdges()]
                if len(incident_indices) != 2:
                    continue
                other_indices = [
                    i for i in incident_indices if i != int(short_edge.index)]
                if len(other_indices) != 1:
                    continue
                other_edge = myPrt.edges[other_indices[0]]
                other_size = float(other_edge.getSize(printResults=False))
                if other_size < threshold:
                    continue
                other_faces = set([int(i) for i in other_edge.getFaces()])
                if other_faces != short_faces:
                    continue

                origin = np.asarray(vertex.pointOn[0], dtype=float)

                def vector_from_vertex(edge):
                    endpoint_indices = [int(i) for i in edge.getVertices()]
                    opposite = [i for i in endpoint_indices if i != vertex_index]
                    if len(opposite) != 1:
                        return None
                    point = np.asarray(
                        myPrt.vertices[opposite[0]].pointOn[0], dtype=float)
                    vector = point - origin
                    length = float(np.linalg.norm(vector))
                    if length <= 0.0:
                        return None
                    return vector / length

                short_vector = vector_from_vertex(short_edge)
                other_vector = vector_from_vertex(other_edge)
                if short_vector is None or other_vector is None:
                    continue
                cosine = float(np.dot(short_vector, other_vector))
                if cosine > cosine_limit:
                    continue

                candidates[vertex_index] = (
                    vertex,
                    int(short_edge.index),
                    float(short_size),
                    int(other_edge.index),
                    float(other_size),
                    cosine,
                )

        if candidates:
            vertices = tuple([
                candidates[index][0] for index in sorted(candidates.keys())])
            myPrt.ignoreEntity(entities=vertices)
            for index in sorted(candidates.keys()):
                record = candidates[index]
                print('>>> SHORT_EDGE_VIRTUAL_TOPOLOGY vertex=%d '
                      'short_edge=%d short_size=%.12g other_edge=%d '
                      'other_size=%.12g cosine=%.12g'
                      % (index, record[1], record[2], record[3],
                         record[4], record[5]))

        remaining_short_edges = []
        for edge in myPrt.edges:
            size = float(edge.getSize(printResults=False))
            if size < threshold:
                remaining_short_edges.append((int(edge.index), size))

        report_path = self.pMesh.get('short_edge_report_path')
        if report_path:
            with open(os.path.abspath(str(report_path)), 'w') as report_file:
                report_file.write('threshold=%.12g\n' % threshold)
                report_file.write(
                    'short_edges_before=%d\n' % len(short_edge_records))
                report_file.write(
                    'ignored_vertices=%d\n' % len(candidates))
                report_file.write(
                    'short_edges_after=%d\n' % len(remaining_short_edges))
                for edge_index, size in short_edge_records:
                    report_file.write(
                        'before_edge=%d size=%.12g\n' % (edge_index, size))
                for index in sorted(candidates.keys()):
                    record = candidates[index]
                    report_file.write(
                        'ignored_vertex=%d short_edge=%d short_size=%.12g '
                        'other_edge=%d other_size=%.12g cosine=%.12g\n'
                        % (index, record[1], record[2], record[3],
                           record[4], record[5]))
                for edge_index, size in remaining_short_edges:
                    report_file.write(
                        'after_edge=%d size=%.12g\n' % (edge_index, size))

    def create_mesh(self):

        myPrt = self.model.parts[self.name_part]
        mesh_control_scope = str(
            self.pMesh.get('shell_mesh_control_scope', 'cells')
        ).lower()
        if mesh_control_scope == 'faces':
            regions = myPrt.faces
            myPrt.setMeshControls(
                regions=regions,
                elemShape=QUAD_DOMINATED,
                technique=FREE,
            )
        elif mesh_control_scope == 'cells':
            regions = myPrt.cells
            myPrt.setMeshControls(
                regions=regions,
                elemShape=QUAD,
                technique=STRUCTURED,
            )
        else:
            raise ValueError(
                'Invalid shell_mesh_control_scope: %s' % mesh_control_scope
            )

        # Some lofted stringer faces are five-sided polygons.  Abaqus' free
        # QUAD_DOMINATED transition can, for isolated station geometries,
        # create two distinct mesh nodes at exactly the same coordinates and
        # hence zero-area S3 transition elements.  Keep the production default
        # unchanged, but provide a controlled pure-triangle strategy for all
        # stringer faces so it can be datachecked before release.
        stringer_control = str(
            self.pMesh.get(
                'stringer_mesh_control', 'free_quad_dominated')
        ).lower()
        if stringer_control in ('free_tri', 'free_quad_medial_axis'):
            controlled_faces = set()
            controlled_sets = 0
            target_sets = self.pMesh.get(
                'stringer_mesh_control_target_sets', None)
            if target_sets is not None:
                target_sets = set([str(value) for value in target_sets])
                missing_sets = sorted([
                    value for value in target_sets
                    if value not in myPrt.sets.keys()
                ])
                if missing_sets:
                    raise ValueError(
                        'Unknown stringer_mesh_control_target_sets: %s'
                        % missing_sets)
            for name_set in sorted(myPrt.sets.keys()):
                if (not str(name_set).startswith('face_wingbox')
                        or '_stringer_' not in str(name_set)):
                    continue
                if target_sets is not None and str(name_set) not in target_sets:
                    continue
                faces = myPrt.sets[name_set].faces
                new_faces = tuple([
                    face for face in faces
                    if int(face.index) not in controlled_faces
                ])
                if not new_faces:
                    continue
                if stringer_control == 'free_tri':
                    myPrt.setMeshControls(
                        regions=new_faces,
                        elemShape=TRI,
                        technique=FREE,
                    )
                else:
                    myPrt.setMeshControls(
                        regions=new_faces,
                        elemShape=QUAD_DOMINATED,
                        technique=FREE,
                        algorithm=MEDIAL_AXIS,
                    )
                for face in new_faces:
                    controlled_faces.add(int(face.index))
                controlled_sets += 1
            print('>>> STRINGER_MESH_CONTROL mode=%s faces=%d sets=%d '
                  'targets=%s'
                  % (stringer_control, len(controlled_faces), controlled_sets,
                     sorted(target_sets) if target_sets is not None else 'all'))
        elif stringer_control != 'free_quad_dominated':
            raise ValueError(
                'Invalid stringer_mesh_control: %s' % stringer_control
            )
        myPrt.generateMesh()

    def set_element_type(self):
        '''
        Set element type as Shell elements (S4R) for all faces.
        '''
        myPrt = self.model.parts[self.name_part]
        elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD, 
            secondOrderAccuracy=OFF, hourglassControl=DEFAULT)
        elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
        myPrt.setElementType(regions=myPrt.sets['all'], elemTypes=(elemType1, elemType2))

    def set_section_assignment(self):
        '''
        Set section assignment via Abaqus Module: Property
        '''
        myPrt = self.model.parts[self.name_part]
        n_sections = len(self.sections)
        
        material_name=str(self.pMesh['material_name'])
        
        for i_section in range(n_sections - 1):
            sec0 = self.sections[i_section]
            tag = 'wingbox%d' % i_section

            # Cover faces
            for side in ['upper', 'lower']:
                params = get_span_group_layup_params(
                    self.pMesh, get_cover_layup_params(self.pMesh, side),
                    i_section)
                
                primaryAxisVector = get_primaryAxisVector_spanwise(
                    self.sections, i_section, feature='cover', side=side)
                
                name_set = 'face_%s_cover_%s' % (tag, side)
                self.name_layups.append(name_set)
                create_shell_CompositeLayup_of_set(
                    myPrt=myPrt, name_set=name_set,
                    ply_thickness=get_design_region_ply_thickness(
                        self.pMesh, params, 'cover', i_section),
                    ply_angles=params['layup_orientAngles'],
                    name_surface=name_set,
                    primaryAxisVector=primaryAxisVector,
                    symmetric=params['layup_symmetric'],
                    numIntPoints=self.pMesh['ply_numIntPts'],
                    material_name=material_name)

            # Spar faces
            params = self.pMesh['spar']
            for j in range(sec0.n_spars):

                params_j = get_span_group_layup_params(
                    self.pMesh, params[j], i_section)

                primaryAxisVector = get_primaryAxisVector_spanwise(
                    self.sections, i_section, feature='spar', index=j)
                
                name_set = 'face_%s_spar_%d' % (tag, j)
                self.name_layups.append(name_set)
                create_shell_CompositeLayup_of_set(
                    myPrt=myPrt, name_set=name_set,
                    ply_thickness=get_design_region_ply_thickness(
                        self.pMesh, params_j, 'spar', i_section),
                    ply_angles=params_j['layup_orientAngles'],
                    name_surface=name_set,
                    primaryAxisVector=primaryAxisVector,
                    symmetric=params_j['layup_symmetric'],
                    numIntPoints=self.pMesh['ply_numIntPts'],
                    material_name=material_name)

            # Stringer faces (web + flange per side, combined into one set)
            params = get_span_group_layup_params(
                self.pMesh, self.pMesh['stringer'], i_section)
            for j in range(sec0.n_stringers):
                for side in ['upper', 'lower']:

                    primaryAxisVector = get_primaryAxisVector_spanwise(
                        self.sections, i_section, feature='stringer', side=side, index=j)
                    
                    name_set='face_%s_stringer_%d_%s' % (tag, j, side)
                    self.name_layups.append(name_set)
                    create_shell_CompositeLayup_of_set(
                        myPrt=myPrt, name_set=name_set,
                        ply_thickness=get_design_region_ply_thickness(
                            self.pMesh, params, 'stringer', i_section),
                        ply_angles=params['layup_orientAngles'],
                        name_surface=name_set,
                        primaryAxisVector=primaryAxisVector,
                        symmetric=params['layup_symmetric'],
                        numIntPoints=self.pMesh['ply_numIntPts'],
                        material_name=material_name)
