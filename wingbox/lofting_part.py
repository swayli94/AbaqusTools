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
        mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=section.chord*10, transform=transform)
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

            transform = myPrt.MakeSketchTransform(
                sketchPlane=self.get_datum_by_name(myPrt, 'XYPLANE-%d'%(i)),
                sketchUpEdge=self.get_datum_by_name(myPrt, 'XAXIS'), 
                sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 0.0))
        
            mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=1000.0, transform=transform)
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
                    edges = self.get_edges(myPrt, findAt_points=findAt_points)
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
                        edges = self.get_edges(myPrt, findAt_points=findAt_points)
                        loftsections.append(tuple(edges))
                    myPrt.ShellLoft(loftsections=tuple(loftsections), startCondition=NONE, endCondition=NONE)
                    self.rename_feature(myPrt, 'wingbox%d_stringer_%d_%s'%(i_section, j, side))
        
        #* Delete the reference plane
        myPrt.setValues(geometryRefinement=EXTRA_FINE)
        del myPrt.features['reference_plane']

    def create_surface(self):
        '''
        Create surfaces, including:
        - cover surfaces (upper/lower)
        - spar surfaces (each spar)
        - stringer surfaces (each stringer, upper/lower)
        '''
        n_sections = len(self.sections)

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
                faces = self.get_faces(myPrt, findAt_points, getClosest=True, searchTolerance=1E-2)
                myPrt.Surface(side1Faces=faces, name='face_%s_cover_%s' % (tag, side))

            # Spar faces
            for j in range(sec0.n_spars):
                pts0 = sec0.get_selection_points(feature='spar', side=None, index=j)
                pts1 = sec1.get_selection_points(feature='spar', side=None, index=j)
                findAt_points=mid_pts(pts0, pts1)
                faces = self.get_faces(myPrt, findAt_points, getClosest=True, searchTolerance=1E-2)
                myPrt.Surface(side1Faces=faces, name='face_%s_spar_%d' % (tag, j))

            # Stringer faces (web + flange per side, combined into one set)
            for j in range(sec0.n_stringers):
                for side in ['upper', 'lower']:
                    pts0 = sec0.get_selection_points(feature='stringer', side=side, index=j)
                    pts1 = sec1.get_selection_points(feature='stringer', side=side, index=j)
                    findAt_points=mid_pts(pts0, pts1)
                    faces = self.get_faces(myPrt, findAt_points, getClosest=True, searchTolerance=1E-2)
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
                    getClosest=True, searchTolerance=1E-2)

            # Spar faces
            for j in range(sec0.n_spars):
                pts0 = sec0.get_selection_points(feature='spar', side=None, index=j)
                pts1 = sec1.get_selection_points(feature='spar', side=None, index=j)
                self.create_geometry_set(
                    name_set='face_%s_spar_%d' % (tag, j),
                    findAt_points=mid_pts(pts0, pts1),
                    geometry='face')

            # Stringer faces (web + flange per side, combined into one set)
            for j in range(sec0.n_stringers):
                for side in ['upper', 'lower']:
                    pts0 = sec0.get_selection_points(feature='stringer', side=side, index=j)
                    pts1 = sec1.get_selection_points(feature='stringer', side=side, index=j)
                    self.create_geometry_set(
                        name_set='face_%s_stringer_%d_%s' % (tag, j, side),
                        findAt_points=mid_pts(pts0, pts1),
                        geometry='face')

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
                    geometry='edge')
            
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

    def create_mesh(self):

        myPrt = self.model.parts[self.name_part]
        myPrt.setMeshControls(regions=myPrt.cells, elemShape=QUAD, technique=STRUCTURED)
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
            params = self.pMesh['cover']
            for side in ['upper', 'lower']:
                
                primaryAxisVector = get_primaryAxisVector_spanwise(
                    self.sections, i_section, feature='cover', side=side)
                
                name_set = 'face_%s_cover_%s' % (tag, side)
                self.name_layups.append(name_set)
                create_shell_CompositeLayup_of_set(
                    myPrt=myPrt, name_set=name_set,
                    ply_thickness=self.pMesh['ply_thickness'],
                    ply_angles=params['layup_orientAngles'],
                    name_surface=name_set,
                    primaryAxisVector=primaryAxisVector,
                    symmetric=params['layup_symmetric'],
                    numIntPoints=self.pMesh['ply_numIntPts'],
                    material_name=material_name)

            # Spar faces
            params = self.pMesh['spar']
            for j in range(sec0.n_spars):

                primaryAxisVector = get_primaryAxisVector_spanwise(
                    self.sections, i_section, feature='spar', index=j)
                
                name_set = 'face_%s_spar_%d' % (tag, j)
                self.name_layups.append(name_set)
                create_shell_CompositeLayup_of_set(
                    myPrt=myPrt, name_set=name_set,
                    ply_thickness=self.pMesh['ply_thickness'],
                    ply_angles=params[j]['layup_orientAngles'],
                    name_surface=name_set,
                    primaryAxisVector=primaryAxisVector,
                    symmetric=params[j]['layup_symmetric'],
                    numIntPoints=self.pMesh['ply_numIntPts'],
                    material_name=material_name)

            # Stringer faces (web + flange per side, combined into one set)
            params = self.pMesh['stringer']
            for j in range(sec0.n_stringers):
                for side in ['upper', 'lower']:

                    primaryAxisVector = get_primaryAxisVector_spanwise(
                        self.sections, i_section, feature='stringer', side=side, index=j)
                    
                    name_set='face_%s_stringer_%d_%s' % (tag, j, side)
                    self.name_layups.append(name_set)
                    create_shell_CompositeLayup_of_set(
                        myPrt=myPrt, name_set=name_set,
                        ply_thickness=self.pMesh['ply_thickness'],
                        ply_angles=params['layup_orientAngles'],
                        name_surface=name_set,
                        primaryAxisVector=primaryAxisVector,
                        symmetric=params['layup_symmetric'],
                        numIntPoints=self.pMesh['ply_numIntPts'],
                        material_name=material_name)

