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

from geometry import WingSectionGeometry


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

    def create_set(self):
        '''
        Create sets for ribs, including:
        - rib faces (each face)
        - rib edges (upper/lower/left/right edge of each rib)
        - cutout edges (upper/lower/left/right edge of each cutout,
        and fillet curve edge of each cutout)

        '''
        section = self.section

        def mid_pt(p0, p1):
            return (0.5*(p0[0]+p1[0]), 0.5*(p0[1]+p1[1]), 0.5*(p0[2]+p1[2]))

        def mid_pts(pts0, pts1):
            return [mid_pt(p0, p1) for p0, p1 in zip(pts0, pts1)]

        # Rib face: mid-height point at the chord midpoint (between spars)
        idx = 1
        rib_face_pt = (
            float(0.5 * (section.x3d_upper_cover[idx] + section.x3d_lower_cover[idx])),
            float(0.5 * (section.y3d_upper_cover[idx] + section.y3d_lower_cover[idx])),
            float(section.zLE))
        self.create_geometry_set(
            name_set='face_rib_%d' % (self.index_rib),
            findAt_points=[rib_face_pt],
            geometry='face',
            getClosest=True, searchTolerance=1E-2)

        # Cover boundary edges (upper and lower splines)
        self.create_geometry_set(
            name_set='edge_rib_%d_cover_upper' % (self.index_rib),
            findAt_points=[(float(section.x3d_upper_cover[idx]),
                            float(section.y3d_upper_cover[idx]),
                            float(section.zLE))],
            geometry='edge')
        self.create_geometry_set(
            name_set='edge_rib_%d_cover_lower' % (self.index_rib),
            findAt_points=[(float(section.x3d_lower_cover[idx]),
                            float(section.y3d_lower_cover[idx]),
                            float(section.zLE))],
            geometry='edge')

        # Front and rear spar boundary edges
        for j, label in [(0, 'front'), (-1, 'rear')]:
            pt = section.spars[j].get_selection_point(feature='spar')
            self.create_geometry_set(
                name_set='edge_rib_%d_spar_%s' % (self.index_rib, label),
                findAt_points=[pt],
                geometry='edge')

        # Cutout straight edges and fillet curves
        for k, cutout in enumerate(section.cutouts):
            for side in ['upper', 'lower', 'left', 'right']:
                pt = cutout.get_selection_point(feature='edge', side=side)
                self.create_geometry_set(
                    name_set='edge_rib_%d_cutout_%d_%s' % (self.index_rib, k, side),
                    findAt_points=[pt],
                    geometry='edge')
            for side in ['upper-left', 'upper-right', 'lower-left', 'lower-right']:
                pt = cutout.get_selection_point(feature='fillet-curve', side=side)
                self.create_geometry_set(
                    name_set='edge_rib_%d_cutout_%d_fillet_%s' % (self.index_rib, k, side.replace('-', '_')),
                    findAt_points=[pt],
                    geometry='edge')

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
