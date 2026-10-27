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

from geometry import WingSectionGeometry, get_primaryAxisVector_section
from layup import create_shell_CompositeLayup_of_set


ALUMINUM_RIB_TYPES = ('aluminum', 'aluminium', 'aluminum_alloy', 'aluminium_alloy')


def get_component_ply_thickness(pMesh, params):
    '''
    Get component-level ply thickness, with old global value as fallback.
    '''
    return params.get('ply_thickness', pMesh['ply_thickness'])


def get_non_design_rib_indices(pMesh):
    '''
    Get rib indices attached to non-design wingbox bays.
    '''
    if 'non_design_rib_indices' in pMesh:
        return set([int(i) for i in pMesh.get('non_design_rib_indices', [])])

    rib_indices = set()
    for i_section in pMesh.get('non_design_wingbox_sections', []):
        i_section = int(i_section)
        rib_indices.add(i_section)
        rib_indices.add(i_section + 1)
    return rib_indices


def get_rib_thickness(pMesh, params, index_rib):
    '''
    Get rib thickness with optional non-design rib override.
    '''
    non_design_region = pMesh.get('non_design_region', {})
    if int(index_rib) in get_non_design_rib_indices(pMesh):
        if 'rib_thickness' in non_design_region:
            return float(non_design_region['rib_thickness'])
    return get_component_ply_thickness(pMesh, params)


def get_rib_material_type(pMesh):
    return str(pMesh.get('rib_material_type', 'composite')).lower()


def get_rib_material_name(pMesh):
    if get_rib_material_type(pMesh) in ALUMINUM_RIB_TYPES:
        return str(pMesh.get('rib_material_name', 'Aluminum-7075'))
    return str(pMesh['material_name'])


def get_rib_partition_mode(pMesh):
    return str(pMesh.get('rib_partition_mode', 'shortest_path')).lower()


def use_rib_face_tie_for_internal_spars(pMesh):
    return get_rib_partition_mode(pMesh) in ('face_tie', 'face-tie')


def get_rib_edge_spar_indices(section, pMesh):
    if use_rib_face_tie_for_internal_spars(pMesh) and section.n_spars > 2:
        return [0, section.n_spars - 1]
    return range(section.n_spars)


class RibPart(Part):
    '''
    Class for wingbox rib part.
    '''
    def __init__(self, model, pGeo, pMesh, index_rib=0):
        super(RibPart,self).__init__(model, pGeo, pMesh)
        self.name_part = 'rib_%d'%(index_rib)
        self.index_rib = index_rib

        section_params = pGeo['sections'][index_rib]
        section_params['airfoil'] = os.path.join(path, section_params['airfoil'])
        self.section = WingSectionGeometry()
        self.section.set_parameters(section_params)
        
        self.name_layups = []

    def print_debug_counts(self, stage):
        '''
        Print compact geometry/mesh counts for diagnosing invalid rib layouts.
        '''
        try:
            myPrt = self.model.parts[self.name_part]
        except Exception:
            print('[rib debug] %s %s: part not available' % (self.name_part, stage))
            return
        print(
            '[rib debug] %s %s: zLE=%.6g faces=%d edges=%d nodes=%d elements=%d'
            % (
                self.name_part,
                stage,
                float(self.section.zLE),
                len(myPrt.faces),
                len(myPrt.edges),
                len(myPrt.nodes),
                len(myPrt.elements),
            )
        )

    def create_sketch(self):
        '''
        Create stand-alone sketches via Abaqus Module: Sketch
        
        - Create Sketch
        '''
        self._create_sketch_rib()
    
    def _create_sketch_rib(self):
        '''
        Create stand-alone sketch for rib: cover + cutout
        '''
        section = self.section
            
        mySkt = self.model.ConstrainedSketch(name=self.name_part, sheetSize=section.chord*3)
        
        cover_u = np.concatenate((section.x3d_upper_cover[:,np.newaxis], section.y3d_upper_cover[:,np.newaxis]), axis=1) # [n,2]
        cover_l = np.concatenate((section.x3d_lower_cover[:,np.newaxis], section.y3d_lower_cover[:,np.newaxis]), axis=1) # [n,2]
        
        mySkt.Spline(points=tuple(cover_u))
        mySkt.Spline(points=tuple(cover_l))
        
        for j in [0, -1]:
            pt1 = section.spars[j].get_selection_point(feature='root', side='upper')[:2]
            pt2 = section.spars[j].get_selection_point(feature='root', side='lower')[:2]
            mySkt.Line(point1=pt1, point2=pt2)

        for j, cutout in enumerate(section.cutouts):
            # edge
            for side1, side2 in [('upper-left', 'upper-right'), ('lower-left', 'lower-right'),
                                    ('left-lower', 'left-upper'), ('right-lower', 'right-upper')]:
                pt1 = cutout.get_selection_point(feature='corner', side=side1)[:2]
                pt2 = cutout.get_selection_point(feature='corner', side=side2)[:2]
                mySkt.Line(point1=pt1, point2=pt2)
            # fillet
            for side_corner, side1, side2 in [
                ('upper-left', 'upper-left', 'left-upper'),
                ('upper-right', 'right-upper', 'upper-right'),
                ('lower-left', 'left-lower', 'lower-left'),
                ('lower-right', 'lower-right', 'right-lower')]:
                pt1 = cutout.get_selection_point(feature='corner', side=side1)[:2]
                pt2 = cutout.get_selection_point(feature='corner', side=side2)[:2]
                pt3 = cutout.get_selection_point(feature='fillet-curve', side=side_corner)[:2]
                mySkt.Arc3Points(point1=pt1, point2=pt2, point3=pt3)
    
    def create_part(self):
        '''
        Create part for rib via Abaqus Module: Part
        '''
        myPrt = self.model.Part(name=self.name_part, dimensionality=THREE_D, type=DEFORMABLE_BODY)

        #* Reference plane and axis
        myPrt.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=self.section.zLE)
        self.rename_feature(myPrt, 'XYPLANE')
        myPrt.DatumAxisByPrincipalAxis(principalAxis=XAXIS)
        self.rename_feature(myPrt, 'XAXIS')
        myPrt.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
        self.rename_feature(myPrt, 'YAXIS')
        myPrt.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
        self.rename_feature(myPrt, 'ZAXIS')
        self.create_datum_csys_3p(myPrt, 'csys_plate', origin=[0.0, 0.0, 0.0],
                                    dx=[1, 0, 0], dy=[0, 1, 0])

        #* Plane for plate sketch
        transform = myPrt.MakeSketchTransform(
            sketchPlane=self.get_datum_by_name(myPrt, 'XYPLANE'),
            sketchUpEdge=self.get_datum_by_name(myPrt, 'XAXIS'), 
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 0.0))
    
        #* Section sketch
        mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=self.section.chord*3, transform=transform)
        mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)
        mySkt.retrieveSketch(sketch=self.model.sketches[self.name_part])
        
        #* Part by Shell (Planer)
        myPrt.BaseShell(sketch=mySkt)

        #* Post procedure
        myPrt.setValues(geometryRefinement=EXTRA_FINE)
        del self.model.sketches['__profile__']
        self.print_debug_counts('after create_part')

    def create_partition(self):
        '''
        Partition rib faces with internal spars.
        '''
        section = self.section
        if section.n_spars <= 2:
            return

        mode = get_rib_partition_mode(self.pMesh)
        if mode in ('none', 'off', 'false', '0', 'face_tie', 'face-tie'):
            return
        if mode == 'shortest_path':
            self._partition_internal_spars_shortest_path()
        elif mode == 'sketch':
            self._partition_internal_spars_by_sketch()
        elif mode in ('datum_plane', 'datum'):
            self._partition_internal_spars_by_datum_plane()
        elif mode == 'hybrid':
            self._partition_internal_spars_hybrid()
        else:
            raise ValueError(
                "Unsupported rib_partition_mode %r for %s. "
                "Use 'shortest_path', 'sketch', 'datum_plane', 'hybrid', 'face_tie', or 'none'."
                % (mode, self.name_part)
            )

    def _validate_partition_faces(self, n_faces_before):
        myPrt = self.model.parts[self.name_part]
        self.print_debug_counts('after create_partition')
        expected_min_faces = n_faces_before + self.section.n_spars - 2
        if len(myPrt.faces) < expected_min_faces:
            raise RuntimeError(
                "Rib partition failed for %s (index=%d, zLE=%.6g): "
                "expected at least %d faces after partition, got %d."
                % (
                    self.name_part,
                    self.index_rib,
                    float(self.section.zLE),
                    expected_min_faces,
                    len(myPrt.faces),
                )
            )

    def _partition_internal_spars_shortest_path(self):
        section = self.section
        myPrt = self.model.parts[self.name_part]
        n_faces_before = len(myPrt.faces)
        for j in range(1, section.n_spars - 1):
            pt1 = section.spars[j].get_selection_point(feature='root', side='upper')
            pt2 = section.spars[j].get_selection_point(feature='root', side='lower')
            myPrt.PartitionFaceByShortestPath(
                faces=myPrt.faces,
                point1=pt1,
                point2=pt2)
        self._validate_partition_faces(n_faces_before)

    def _partition_internal_spars_by_sketch(self):
        section = self.section
        myPrt = self.model.parts[self.name_part]
        n_faces_before = len(myPrt.faces)

        transform = myPrt.MakeSketchTransform(
            sketchPlane=self.get_datum_by_name(myPrt, 'XYPLANE'),
            sketchUpEdge=self.get_datum_by_name(myPrt, 'XAXIS'),
            sketchPlaneSide=SIDE1,
            sketchOrientation=BOTTOM,
            origin=(0.0, 0.0, 0.0))
        mySkt = self.model.ConstrainedSketch(
            name='__rib_partition__',
            sheetSize=section.chord * 3,
            transform=transform)
        mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)
        for j in range(1, section.n_spars - 1):
            pt1 = section.spars[j].get_selection_point(feature='root', side='upper')[:2]
            pt2 = section.spars[j].get_selection_point(feature='root', side='lower')[:2]
            mySkt.Line(point1=pt1, point2=pt2)

        try:
            myPrt.PartitionFaceBySketch(faces=myPrt.faces, sketch=mySkt)
        finally:
            if '__rib_partition__' in self.model.sketches:
                del self.model.sketches['__rib_partition__']

        self._validate_partition_faces(n_faces_before)

    def _partition_internal_spars_by_datum_plane(self):
        section = self.section
        myPrt = self.model.parts[self.name_part]
        n_faces_before = len(myPrt.faces)

        for j in range(1, section.n_spars - 1):
            pt1 = section.spars[j].get_selection_point(feature='root', side='upper')
            pt2 = section.spars[j].get_selection_point(feature='root', side='lower')
            pt3 = (float(pt1[0]), float(pt1[1]), float(pt1[2]) + 1.0)
            feature = myPrt.DatumPlaneByThreePoints(point1=pt1, point2=pt2, point3=pt3)
            datum = myPrt.datums[feature.id]
            myPrt.PartitionFaceByDatumPlane(faces=myPrt.faces, datumPlane=datum)

        self._validate_partition_faces(n_faces_before)

    def _partition_internal_spars_hybrid(self):
        errors = []
        for mode_name, partition_method in [
                ('datum_plane', self._partition_internal_spars_by_datum_plane),
                ('sketch', self._partition_internal_spars_by_sketch),
                ('shortest_path', self._partition_internal_spars_shortest_path)]:
            try:
                partition_method()
                return
            except Exception as exc:
                errors.append('%s: %s' % (mode_name, exc))

        raise RuntimeError(
            "Rib hybrid partition failed for %s (index=%d, zLE=%.6g). Attempts: %s"
            % (
                self.name_part,
                self.index_rib,
                float(self.section.zLE),
                " | ".join(errors),
            )
        )

    def create_surface(self):
        '''
        Create surfaces
        '''
        myPrt = self.model.parts[self.name_part]
        myPrt.Surface(side1Faces=myPrt.faces, name='face_rib')

    def create_set(self):
        '''
        Create sets for ribs, including:
        - rib faces (each face)
        - rib edges (upper/lower/left/right edge of each rib)
        - cutout edges (upper/lower/left/right edge of each cutout,
        and fillet curve edge of each cutout)

        '''
        section = self.section

        myPrt = self.model.parts[self.name_part]
        myPrt.Set(faces=myPrt.faces, name='all') 

        myPrt.Set(faces=myPrt.faces, name='face_rib')

        # Cover boundary edges (upper and lower splines)
        idx = 1
        self.create_geometry_set(
            name_set='edge_rib_cover_upper',
            findAt_points=[(float(section.x3d_upper_cover[idx]),
                            float(section.y3d_upper_cover[idx]),
                            float(section.zLE))],
            geometry='edge')
        self.create_geometry_set(
            name_set='edge_rib_cover_lower',
            findAt_points=[(float(section.x3d_lower_cover[idx]),
                            float(section.y3d_lower_cover[idx]),
                            float(section.zLE))],
            geometry='edge')

        # Front and rear spar boundary edges
        for j in get_rib_edge_spar_indices(self.section, self.pMesh):
            pt = section.spars[j].get_selection_point(feature='spar')
            self.create_geometry_set(
                name_set='edge_rib_spar_%d' % (j),
                findAt_points=[pt],
                geometry='edge')

        # Cutout straight edges and fillet curves
        for k, cutout in enumerate(section.cutouts):
            for side in ['upper', 'lower', 'left', 'right']:
                pt = cutout.get_selection_point(feature='edge', side=side)
                self.create_geometry_set(
                    name_set='edge_rib_cutout_%d_%s' % (k, side),
                    findAt_points=[pt],
                    geometry='edge')
            for side in ['upper-left', 'upper-right', 'lower-left', 'lower-right']:
                pt = cutout.get_selection_point(feature='fillet-curve', side=side)
                self.create_geometry_set(
                    name_set='edge_rib_cutout_%d_fillet_%s' % (k, side.replace('-', '_')),
                    findAt_points=[pt],
                    geometry='edge')

        self.validate_required_sets()

    def validate_required_sets(self):
        '''
        Check required rib sets immediately after set creation.

        Missing rib sets otherwise appear much later as Abaqus input-file
        surface errors, which makes invalid rib layouts hard to diagnose.
        '''
        myPrt = self.model.parts[self.name_part]
        required = [('face_rib', 'faces')]
        required.extend([
            ('edge_rib_cover_upper', 'edges'),
            ('edge_rib_cover_lower', 'edges'),
        ])
        for j in get_rib_edge_spar_indices(self.section, self.pMesh):
            required.append(('edge_rib_spar_%d' % j, 'edges'))

        missing = []
        empty = []
        for name_set, attribute in required:
            if name_set not in myPrt.sets:
                missing.append(name_set)
                continue
            geometry_array = getattr(myPrt.sets[name_set], attribute)
            if len(geometry_array) == 0:
                empty.append(name_set)

        if missing or empty:
            raise RuntimeError(
                "Rib set generation failed for %s (index=%d, zLE=%.6g). "
                "Missing sets: %s. Empty sets: %s."
                % (
                    self.name_part,
                    self.index_rib,
                    float(self.section.zLE),
                    missing,
                    empty,
                )
            )

    def validate_required_mesh_sets(self):
        '''
        Check required rib sets after mesh generation.

        Abaqus can accept geometry sets in CAE but later fail while writing the
        input deck if those sets do not resolve to mesh nodes/elements.
        '''
        myPrt = self.model.parts[self.name_part]
        if len(myPrt.nodes) == 0 or len(myPrt.elements) == 0:
            raise RuntimeError(
                "Rib mesh generation failed for %s (index=%d, zLE=%.6g): "
                "nodes=%d, elements=%d."
                % (
                    self.name_part,
                    self.index_rib,
                    float(self.section.zLE),
                    len(myPrt.nodes),
                    len(myPrt.elements),
                )
            )

        required = [('face_rib', 'elements')]
        required.extend([
            ('edge_rib_cover_upper', 'nodes'),
            ('edge_rib_cover_lower', 'nodes'),
        ])
        for j in get_rib_edge_spar_indices(self.section, self.pMesh):
            required.append(('edge_rib_spar_%d' % j, 'nodes'))

        missing = []
        empty = []
        unavailable = []
        for name_set, attribute in required:
            if name_set not in myPrt.sets:
                missing.append(name_set)
                continue
            try:
                mesh_array = getattr(myPrt.sets[name_set], attribute)
                count = len(mesh_array)
            except Exception:
                unavailable.append('%s.%s' % (name_set, attribute))
                continue
            if count == 0:
                empty.append('%s.%s' % (name_set, attribute))

        if missing or empty or unavailable:
            raise RuntimeError(
                "Rib mesh set validation failed for %s (index=%d, zLE=%.6g). "
                "Missing sets: %s. Empty mesh sets: %s. Unavailable mesh arrays: %s."
                % (
                    self.name_part,
                    self.index_rib,
                    float(self.section.zLE),
                    missing,
                    empty,
                    unavailable,
                )
            )

    def set_seeding(self):
        
        myPrt = self.model.parts[self.name_part]
        myPrt.seedPart(size=self.pMesh['rib_seedPart_size'], 
                        deviationFactor=0.1, minSizeFactor=0.1)

        # Spar side edges (for seeding)
        for j in get_rib_edge_spar_indices(self.section, self.pMesh):
            name_set='edge_rib_spar_%d' % (j)
            myPrt.seedEdgeBySize(edges=myPrt.sets[name_set].edges,
                    size=self.pMesh['rib_seedEdge_size'],
                    deviationFactor=0.1, constraint=FINER)
                
        # Cutout edges (for seeding)
        for k, cutout in enumerate(self.section.cutouts):
            for side in ['upper', 'lower', 'left', 'right']:
                name_set = 'edge_rib_cutout_%d_%s' % (k, side)
                myPrt.seedEdgeBySize(edges=myPrt.sets[name_set].edges,
                    size=self.pMesh['rib_seedEdge_size'],
                    deviationFactor=0.1, constraint=FINER)
            for side in ['upper-left', 'upper-right', 'lower-left', 'lower-right']:
                name_set='edge_rib_cutout_%d_fillet_%s' % (k, side.replace('-', '_'))
                myPrt.seedEdgeBySize(edges=myPrt.sets[name_set].edges,
                    size=self.pMesh['fillet_seedEdge_size'],
                    deviationFactor=0.1, constraint=FINER)

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
            myPrt.setMeshControls(regions=regions, elemShape=QUAD_DOMINATED)
        else:
            raise ValueError(
                'Invalid shell_mesh_control_scope: %s' % mesh_control_scope
            )
        myPrt.generateMesh()
        self.print_debug_counts('after create_mesh')
        self.validate_required_mesh_sets()
    
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
        params = self.pMesh['rib'][self.index_rib]
        material_name = get_rib_material_name(self.pMesh)

        if get_rib_material_type(self.pMesh) in ALUMINUM_RIB_TYPES:
            section_name = 'section_%s' % self.name_part
            self.model.HomogeneousShellSection(
                name=section_name,
                preIntegrate=OFF,
                material=material_name,
                thicknessType=UNIFORM,
                thickness=get_rib_thickness(self.pMesh, params, self.index_rib),
                thicknessField='',
                nodalThicknessField='',
                idealization=NO_IDEALIZATION,
                poissonDefinition=DEFAULT,
                thicknessModulus=None,
                temperature=GRADIENT,
                useDensity=OFF,
                integrationRule=SIMPSON,
                numIntPts=self.pMesh.get('rib_numIntPts', 5))
            myPrt.SectionAssignment(
                region=myPrt.sets['face_rib'],
                sectionName=section_name,
                offset=0.0,
                offsetType=MIDDLE_SURFACE,
                offsetField='',
                thicknessAssignment=FROM_SECTION)
            return

        primaryAxisVector = get_primaryAxisVector_section(
            self.section, feature='rib')
        
        name_set = 'face_rib'
        self.name_layups.append(name_set)
        create_shell_CompositeLayup_of_set(
            myPrt=myPrt, name_set=name_set,
            ply_thickness=get_rib_thickness(self.pMesh, params, self.index_rib),
            ply_angles=params['layup_orientAngles'],
            name_surface=name_set,
            primaryAxisVector=primaryAxisVector,
            symmetric=params['layup_symmetric'],
            numIntPoints=self.pMesh['ply_numIntPts'],
            material_name=material_name)
        
