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
            
        mySkt = self.model.ConstrainedSketch(name=self.name_part, sheetSize=1000.0)
        
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
        mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=1000.0, transform=transform)
        mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)
        mySkt.retrieveSketch(sketch=self.model.sketches[self.name_part])
        
        #* Part by Shell (Planer)
        myPrt.BaseShell(sketch=mySkt)

        #* Post procedure
        myPrt.setValues(geometryRefinement=EXTRA_FINE)
        del self.model.sketches['__profile__']

    def create_surface(self):
        '''
        Create surfaces
        '''
        myPrt = self.model.parts[self.name_part]
        section = self.section
        idx = 1
        rib_face_pt = (
            float(0.5 * (section.x3d_upper_cover[idx] + section.x3d_lower_cover[idx])),
            float(0.5 * (section.y3d_upper_cover[idx] + section.y3d_lower_cover[idx])),
            float(section.zLE))
        faces = self.get_faces(myPrt, [rib_face_pt], getClosest=True, searchTolerance=1E-2)
        myPrt.Surface(side1Faces=faces, name='face_rib')

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

        # Rib face: mid-height point at the chord midpoint (between spars)
        idx = 1
        rib_face_pt = (
            float(0.5 * (section.x3d_upper_cover[idx] + section.x3d_lower_cover[idx])),
            float(0.5 * (section.y3d_upper_cover[idx] + section.y3d_lower_cover[idx])),
            float(section.zLE))
        self.create_geometry_set(
            name_set='face_rib',
            findAt_points=[rib_face_pt],
            geometry='face',
            getClosest=True, searchTolerance=1E-2)

        # Cover boundary edges (upper and lower splines)
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
        for j in range(self.section.n_spars):
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

    def set_seeding(self):
        
        myPrt = self.model.parts[self.name_part]
        myPrt.seedPart(size=self.pMesh['rib_seedPart_size'], 
                        deviationFactor=0.1, minSizeFactor=0.1)

        # Spar side edges (for seeding)
        for j in range(self.section.n_spars):
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
        myPrt.setMeshControls(regions=myPrt.cells, elemShape=QUAD_DOMINATED)
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
        params = self.pMesh['rib'][self.index_rib]
        material_name=str(self.pMesh['material_name'])

        primaryAxisVector = get_primaryAxisVector_section(
            self.section, feature='rib')
        
        name_set = 'face_rib'
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
        