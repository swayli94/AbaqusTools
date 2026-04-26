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


class WingboxPart(Part):
    '''
    Class for wingbox part.
    '''
    def __init__(self, name_part, model, pGeo, pMesh):
        super(WingboxPart,self).__init__(model, pGeo, pMesh)
        self.name_part = name_part

        self.sections = []
        for section_params in pGeo['sections']:
            section_params['airfoil'] = os.path.join(path, section_params['airfoil'])
            wsg = WingSectionGeometry()
            wsg.set_parameters(section_params)
            self.sections.append(wsg)
            print('WingboxPart: section %d initialized.'%(len(self.sections)))

    def create_sketch(self):
        '''
        Create stand-alone sketches via Abaqus Module: Sketch
        
        - Create Sketch
        '''
        self._create_sketch_rib()
        self._create_sketch_section_for_sweep()
        self._create_sketch_reference_plane()
    
    def _create_sketch_rib(self):
        '''
        Create stand-alone sketch for rib: cover + cutout
        '''
        #* Create the open-hole plate sketch in x-y plane
        for i, section in enumerate(self.sections):
            
            mySkt = self.model.ConstrainedSketch(name='rib_%d'%(i), sheetSize=section.chord*3)
            
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
        
        for i in range(len(self.sections)):
            self._create_part_rib(i)
            
        self._create_part_lofting()

    def _create_part_rib(self, i_section):
        '''
        Create part for rib via Abaqus Module: Part
        '''
        section = self.sections[i_section]
        name = 'rib_%d'%(i_section)
        
        myPrt = self.model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

        #* Reference plane and axis
        myPrt.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=section.zLE)
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
        mySkt = self.model.ConstrainedSketch(name='__profile__', sheetSize=section.chord*10, transform=transform)
        mySkt.sketchOptions.setValues(gridOrigin=(0.0, 0.0), gridAngle=0.0)
        mySkt.retrieveSketch(sketch=self.model.sketches[name])
        
        #* Part by Shell (Planer)
        myPrt.BaseShell(sketch=mySkt)

        #* Post procedure
        myPrt.setValues(geometryRefinement=EXTRA_FINE)
        del self.model.sketches['__profile__']

    def _create_part_lofting(self):
        '''
        Create part by shell lofting
        '''
        myPrt = self.model.Part(name='lofted_part', dimensionality=THREE_D, type=DEFORMABLE_BODY)

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
