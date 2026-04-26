'''
Geometry description of the wing box structures
'''
import os

import numpy as np
from scipy.interpolate import interp1d

from utils import transform_curve, read_airfoil, split_upper_lower, reconstruct_airfoil


class SparGeometry:
    '''
    Geometry description of the spar.
    
    Spar is a vertical structural member in the wing box,
    which is modelled as a line in the x-y plane.
    
    Parameters
    ---------------
    x: float
        x coordinate of the spar (relative to chord)
    t: float
        thickness of the spar (mm)
    
    Attributes
    ---------------
    n_points: int
        number of points to discretize the spar line.
    '''
    def __init__(self, x, t):
        self._x = x
        self._t = t
        self._n_points = 11
        self._i_mid_point = self._n_points // 2
        
    def update_attributes(self, y_upper, y_lower,
                        chord, twist_degrees, xLE, yLE, zLE):
        '''
        Update the attributes of the spar geometry.
        
        Parameters
        ---------------
        y_upper, y_lower: float
            y coordinates of the upper and lower surfaces of the airfoil
            at the x coordinate of the spar (relative to chord)
        '''
        self.xx = np.linspace(self._x, self._x, num=self._n_points, endpoint=True)
        self.yy = np.linspace(y_lower, y_upper, num=self._n_points, endpoint=True)
        
        self.x3d, self.y3d = transform_curve(self.xx, self.yy,
                scale=chord, rotation=twist_degrees, dx=xLE, dy=yLE, center=[0,0])
        self.z3d = np.ones(self._n_points) * zLE
    
    @property
    def x(self):
        '''
        x coordinate of the spar (relative to chord)
        '''
        return float(self._x)

    @property
    def t(self):
        '''
        thickness of the spar (mm)
        '''
        return float(self._t)

    def get_selection_point(self, feature='spar', side='upper'):
        '''
        Get the selection point of features of the stringer.
        
        Parameters
        ---------------
        feature: str
            'spar' or 'root', indicating the feature of the stringer.
        side: str
            'upper' or 'lower', indicating the side of the stringer.
        
        Returns
        ---------------
        point: tuple
            Tuple of (x, y, z) coordinates (mm) of the mid point of the feature line.
        '''
        if feature == 'spar':
            return (self.x3d[self._i_mid_point],
                    self.y3d[self._i_mid_point],
                    self.z3d[self._i_mid_point])
        elif feature == 'root' and side == 'upper':
            return (self.x3d[-1], self.y3d[-1], self.z3d[-1])
        elif feature == 'root' and side == 'lower':
            return (self.x3d[0], self.y3d[0], self.z3d[0])
        else:
            raise ValueError('Invalid side or feature.', side, feature)


class StringerGeometry:
    '''
    Geometry description of the stringer.
    
    Stringer is a structural member attached to the wing cover,
    which is modelled as a T/C-shaped section in the x-y plane.
    The stringers are paired on the upper and lower wing covers.
    
    T-shaped stringer is modelled as a vertical line (web).
    C-shaped stringer is modelled as a vertical line (web) with a horizontal line (flange) on top.
    
    Parameters
    ---------------
    x: float
        x coordinate of the stringer (relative to chord)
    t: float
        thickness of the stringer (mm)
    h: float
        height of the vertical section (web) of the stringer (mm)
    w: float
        width of the horizontal section (flange) of the stringer (mm)
    '''
    def __init__(self, x, t, h, w):
        self._x = x
        self._t = t
        self._h = h
        self._w = w
    
        self._n_points = 11
        self._i_mid_point = self._n_points // 2
            
    def update_attributes(self, y_upper, y_lower,
                        chord, twist_degrees, xLE, yLE, zLE):
        '''
        Update the attributes of the stringer geometry.
        
        Parameters
        ---------------
        y_upper, y_lower: float
            y coordinates of the upper and lower surfaces of the airfoil
            at the x coordinate of the spar (relative to chord)
        '''
        self._update_web_attributes(y_upper, y_lower, chord, twist_degrees, xLE, yLE, zLE)
        self._update_flange_attributes(y_upper, y_lower, chord, twist_degrees, xLE, yLE, zLE)
    
    def _update_web_attributes(self, y_upper, y_lower,
                            chord, twist_degrees, xLE, yLE, zLE):
        '''
        Update the attributes of the web geometry (from root to corner).
        '''
        self.xx_upper_web = np.linspace(self._x, self._x, num=self._n_points, endpoint=True)
        self.yy_upper_web = np.linspace(y_upper, y_upper-self._h, num=self._n_points, endpoint=True)
        
        self.xx_lower_web = np.linspace(self._x, self._x, num=self._n_points, endpoint=True)
        self.yy_lower_web = np.linspace(y_lower, y_lower+self._h, num=self._n_points, endpoint=True)
    
        self.x3d_upper_web, self.y3d_upper_web = transform_curve(self.xx_upper_web, self.yy_upper_web,
                scale=chord, rotation=twist_degrees, dx=xLE, dy=yLE, center=[0,0])
        self.z3d_upper_web = np.ones(self._n_points) * zLE
        
        self.x3d_lower_web, self.y3d_lower_web = transform_curve(self.xx_lower_web, self.yy_lower_web,
                scale=chord, rotation=twist_degrees, dx=xLE, dy=yLE, center=[0,0])
        self.z3d_lower_web = np.ones(self._n_points) * zLE
        
    def _update_flange_attributes(self, y_upper, y_lower,
                            chord, twist_degrees, xLE, yLE, zLE):
        '''
        Update the attributes of the flange geometry (in the increasing x-direction).
        '''
        self.xx_upper_flange = np.linspace(self._x, self._x+self._w, num=self._n_points, endpoint=True)
        self.yy_upper_flange = np.ones(self._n_points) * (y_upper - self._h)
        
        self.xx_lower_flange = np.linspace(self._x, self._x+self._w, num=self._n_points, endpoint=True)
        self.yy_lower_flange = np.ones(self._n_points) * (y_lower + self._h)
                
        self.x3d_upper_flange, self.y3d_upper_flange = transform_curve(self.xx_upper_flange, self.yy_upper_flange,
                scale=chord, rotation=twist_degrees, dx=xLE, dy=yLE, center=[0,0])
        self.z3d_upper_flange = np.ones(self._n_points) * zLE
        
        self.x3d_lower_flange, self.y3d_lower_flange = transform_curve(self.xx_lower_flange, self.yy_lower_flange,
                scale=chord, rotation=twist_degrees, dx=xLE, dy=yLE, center=[0,0])
        self.z3d_lower_flange = np.ones(self._n_points) * zLE
    
    @property
    def x(self):
        '''
        x coordinate of the stringer (relative to chord)
        '''
        return float(self._x)

    @property
    def t(self):
        '''
        thickness of the stringer (mm)
        '''
        return float(self._t)

    @property
    def h(self):
        '''
        height of the vertical section of the stringer (mm)
        '''
        return float(self._h)
    
    @property
    def w(self):
        '''
        width of the horizontal section of the stringer (mm)
        '''
        return float(self._w)

    def get_selection_point(self, feature='web', side='upper'):
        '''
        Get the selection point of features of the stringer.
        
        Parameters
        ---------------
        feature: str
            'web', 'flange', 'root', 'corner' or 'tip',
            indicating the feature of the stringer.
        side: str
            'upper' or 'lower', indicating the side of the stringer.
        
        Returns
        ---------------
        point: tuple
            Tuple of (x, y, z) coordinates (mm) of the mid point of the feature line.
        '''
        if side == 'upper' and feature == 'web':
            return (self.x3d_upper_web[self._i_mid_point],
                    self.y3d_upper_web[self._i_mid_point],
                    self.z3d_upper_web[self._i_mid_point])
        elif side == 'upper' and feature == 'flange':
            return (self.x3d_upper_flange[self._i_mid_point],
                    self.y3d_upper_flange[self._i_mid_point],
                    self.z3d_upper_flange[self._i_mid_point])
        elif side == 'upper' and feature == 'root':
            return (self.x3d_upper_web[0], self.y3d_upper_web[0], self.z3d_upper_web[0])
        elif side == 'upper' and feature == 'corner':
            return (self.x3d_upper_web[-1], self.y3d_upper_web[-1], self.z3d_upper_web[-1])
        elif side == 'upper' and feature == 'tip':
            return (self.x3d_upper_flange[-1], self.y3d_upper_flange[-1], self.z3d_upper_flange[-1])
            
        elif side == 'lower' and feature == 'web':
            return (self.x3d_lower_web[self._i_mid_point],
                    self.y3d_lower_web[self._i_mid_point],
                    self.z3d_lower_web[self._i_mid_point])
        elif side == 'lower' and feature == 'flange':
            return (self.x3d_lower_flange[self._i_mid_point],
                    self.y3d_lower_flange[self._i_mid_point],
                    self.z3d_lower_flange[self._i_mid_point])
        elif side == 'lower' and feature == 'root':
            return (self.x3d_lower_web[0], self.y3d_lower_web[0], self.z3d_lower_web[0])
        elif side == 'lower' and feature == 'corner':
            return (self.x3d_lower_web[-1], self.y3d_lower_web[-1], self.z3d_lower_web[-1])
        elif side == 'lower' and feature == 'tip':
            return (self.x3d_lower_flange[-1], self.y3d_lower_flange[-1], self.z3d_lower_flange[-1])
        
        else:
            raise ValueError('Invalid side or feature.', side, feature)


class CutoutGeometry:
    '''
    Geometry description of the cutout.
    
    Cutout is a hole in the wing rib,
    which is modelled as a rectangle with filleted corners in the x-y plane.
    
    Parameters
    ---------------
    x: float
        x coordinate of the cutout center (relative to chord)
    y: float
        y coordinate of the cutout center (relative to chord)
    lx: float
        length of the cutout in the x-direction (mm)
    ly: float
        length of the cutout in the y-direction (mm)
    r_fillet: float
        radius of the fillet at the cutout corners (mm)
    '''
    def __init__(self, x, y, lx, ly, r_fillet):
        self._x = x
        self._y = y
        self._lx = lx
        self._ly = ly
        self._r_fillet = r_fillet
    
        self._n_points = 11
        self._i_mid_point = self._n_points // 2
    
    def update_attributes(self, y_upper, y_lower,
                        chord, twist_degrees, xLE, yLE, zLE):
        '''
        Update the attributes of the cutout geometry.
        '''
        hlx = self._lx * 0.5
        hly = self._ly * 0.5
        r   = self._r_fillet
        cx  = self._x
        cy  = self._y
        n   = self._n_points

        # Four straight edges (excluding fillet regions)
        self.xx_upper_edge = np.linspace(cx - hlx + r, cx + hlx - r, num=n)
        self.yy_upper_edge = np.ones(n) * (cy + hly)

        self.xx_lower_edge = np.linspace(cx - hlx + r, cx + hlx - r, num=n)
        self.yy_lower_edge = np.ones(n) * (cy - hly)

        self.xx_left_edge  = np.ones(n) * (cx - hlx)
        self.yy_left_edge  = np.linspace(cy - hly + r, cy + hly - r, num=n)

        self.xx_right_edge = np.ones(n) * (cx + hlx)
        self.yy_right_edge = np.linspace(cy - hly + r, cy + hly - r, num=n)

        # Four fillet quarter-circle curves
        theta_ur = np.linspace(0,            np.pi / 2,     num=n)
        theta_ul = np.linspace(np.pi / 2,    np.pi,         num=n)
        theta_ll = np.linspace(np.pi,        3 * np.pi / 2, num=n)
        theta_lr = np.linspace(3 * np.pi / 2, 2 * np.pi,   num=n)

        self.xx_upper_right_fillet = (cx + hlx - r) + r * np.cos(theta_ur)
        self.yy_upper_right_fillet = (cy + hly - r) + r * np.sin(theta_ur)

        self.xx_upper_left_fillet  = (cx - hlx + r) + r * np.cos(theta_ul)
        self.yy_upper_left_fillet  = (cy + hly - r) + r * np.sin(theta_ul)

        self.xx_lower_left_fillet  = (cx - hlx + r) + r * np.cos(theta_ll)
        self.yy_lower_left_fillet  = (cy - hly + r) + r * np.sin(theta_ll)

        self.xx_lower_right_fillet = (cx + hlx - r) + r * np.cos(theta_lr)
        self.yy_lower_right_fillet = (cy - hly + r) + r * np.sin(theta_lr)

        # Four fillet centers (single points stored as length-1 arrays)
        self.xx_upper_left_center  = np.array([cx - hlx + r])
        self.yy_upper_left_center  = np.array([cy + hly - r])

        self.xx_upper_right_center = np.array([cx + hlx - r])
        self.yy_upper_right_center = np.array([cy + hly - r])

        self.xx_lower_left_center  = np.array([cx - hlx + r])
        self.yy_lower_left_center  = np.array([cy - hly + r])

        self.xx_lower_right_center = np.array([cx + hlx - r])
        self.yy_lower_right_center = np.array([cy - hly + r])

        # Transform all curves to 3D
        def _tc(xx, yy):
            return transform_curve(xx, yy,
                    scale=chord, rotation=twist_degrees, dx=xLE, dy=yLE, center=[0, 0])

        self.x3d_upper_edge, self.y3d_upper_edge = _tc(self.xx_upper_edge, self.yy_upper_edge)
        self.z3d_upper_edge = np.ones(n) * zLE

        self.x3d_lower_edge, self.y3d_lower_edge = _tc(self.xx_lower_edge, self.yy_lower_edge)
        self.z3d_lower_edge = np.ones(n) * zLE

        self.x3d_left_edge, self.y3d_left_edge = _tc(self.xx_left_edge, self.yy_left_edge)
        self.z3d_left_edge = np.ones(n) * zLE

        self.x3d_right_edge, self.y3d_right_edge = _tc(self.xx_right_edge, self.yy_right_edge)
        self.z3d_right_edge = np.ones(n) * zLE

        self.x3d_upper_right_fillet, self.y3d_upper_right_fillet = _tc(self.xx_upper_right_fillet, self.yy_upper_right_fillet)
        self.z3d_upper_right_fillet = np.ones(n) * zLE

        self.x3d_upper_left_fillet, self.y3d_upper_left_fillet = _tc(self.xx_upper_left_fillet, self.yy_upper_left_fillet)
        self.z3d_upper_left_fillet = np.ones(n) * zLE

        self.x3d_lower_left_fillet, self.y3d_lower_left_fillet = _tc(self.xx_lower_left_fillet, self.yy_lower_left_fillet)
        self.z3d_lower_left_fillet = np.ones(n) * zLE

        self.x3d_lower_right_fillet, self.y3d_lower_right_fillet = _tc(self.xx_lower_right_fillet, self.yy_lower_right_fillet)
        self.z3d_lower_right_fillet = np.ones(n) * zLE

        self.x3d_upper_left_center, self.y3d_upper_left_center = _tc(self.xx_upper_left_center, self.yy_upper_left_center)
        self.z3d_upper_left_center = np.ones(1) * zLE

        self.x3d_upper_right_center, self.y3d_upper_right_center = _tc(self.xx_upper_right_center, self.yy_upper_right_center)
        self.z3d_upper_right_center = np.ones(1) * zLE

        self.x3d_lower_left_center, self.y3d_lower_left_center = _tc(self.xx_lower_left_center, self.yy_lower_left_center)
        self.z3d_lower_left_center = np.ones(1) * zLE

        self.x3d_lower_right_center, self.y3d_lower_right_center = _tc(self.xx_lower_right_center, self.yy_lower_right_center)
        self.z3d_lower_right_center = np.ones(1) * zLE

    @property
    def x(self):
        '''
        x coordinate of the cutout center (relative to chord)
        '''
        return float(self._x)

    @property
    def y(self):
        '''
        y coordinate of the cutout center (relative to chord)
        '''
        return float(self._y)

    @property
    def lx(self):
        '''
        length of the cutout in the x-direction (mm)
        '''
        return float(self._lx)
    
    @property
    def ly(self):
        '''
        length of the cutout in the y-direction (mm)
        '''
        return float(self._ly)
    
    @property
    def r_fillet(self):
        '''
        radius of the fillet at the cutout corners (mm)
        '''
        return float(self._r_fillet)

    def get_selection_point(self, feature='web', side='upper'):
        '''
        Get the selection point of features of the cutout.
        
        Parameters
        ---------------
        feature: str
            'edge', 'fillet-curve' or 'fillet-center',
            indicating the feature of the cutout.
        side: str
            'upper', 'lower', 'left', 'right',
            'upper-left', 'upper-right', 'lower-left' or 'lower-right',
            indicating the side of the feature.
        
        Returns
        ---------------
        point: tuple
            Tuple of (x, y, z) coordinates (mm) of the mid point of the feature line.
        '''
        i = self._i_mid_point
        if feature == 'edge':
            if side == 'upper':
                return (self.x3d_upper_edge[i], self.y3d_upper_edge[i], self.z3d_upper_edge[i])
            elif side == 'lower':
                return (self.x3d_lower_edge[i], self.y3d_lower_edge[i], self.z3d_lower_edge[i])
            elif side == 'left':
                return (self.x3d_left_edge[i], self.y3d_left_edge[i], self.z3d_left_edge[i])
            elif side == 'right':
                return (self.x3d_right_edge[i], self.y3d_right_edge[i], self.z3d_right_edge[i])
        elif feature == 'corner':
            if side == 'upper-left':
                return (self.x3d_upper_edge[0], self.y3d_upper_edge[0], self.z3d_upper_edge[0])
            elif side == 'upper-right':
                return (self.x3d_upper_edge[-1], self.y3d_upper_edge[-1], self.z3d_upper_edge[-1])
            elif side == 'lower-left':
                return (self.x3d_lower_edge[0], self.y3d_lower_edge[0], self.z3d_lower_edge[0])
            elif side == 'lower-right':
                return (self.x3d_lower_edge[-1], self.y3d_lower_edge[-1], self.z3d_lower_edge[-1])
            elif side == 'left-upper':
                return (self.x3d_left_edge[-1], self.y3d_left_edge[-1], self.z3d_left_edge[-1])
            elif side == 'left-lower':
                return (self.x3d_left_edge[0], self.y3d_left_edge[0], self.z3d_left_edge[0])
            elif side == 'right-upper':
                return (self.x3d_right_edge[-1], self.y3d_right_edge[-1], self.z3d_right_edge[-1])
            elif side == 'right-lower':
                return (self.x3d_right_edge[0], self.y3d_right_edge[0], self.z3d_right_edge[0])
        elif feature == 'fillet-curve':
            if side == 'upper-left':
                return (self.x3d_upper_left_fillet[i], self.y3d_upper_left_fillet[i], self.z3d_upper_left_fillet[i])
            elif side == 'upper-right':
                return (self.x3d_upper_right_fillet[i], self.y3d_upper_right_fillet[i], self.z3d_upper_right_fillet[i])
            elif side == 'lower-left':
                return (self.x3d_lower_left_fillet[i], self.y3d_lower_left_fillet[i], self.z3d_lower_left_fillet[i])
            elif side == 'lower-right':
                return (self.x3d_lower_right_fillet[i], self.y3d_lower_right_fillet[i], self.z3d_lower_right_fillet[i])
        elif feature == 'fillet-center':
            if side == 'upper-left':
                return (self.x3d_upper_left_center[0], self.y3d_upper_left_center[0], self.z3d_upper_left_center[0])
            elif side == 'upper-right':
                return (self.x3d_upper_right_center[0], self.y3d_upper_right_center[0], self.z3d_upper_right_center[0])
            elif side == 'lower-left':
                return (self.x3d_lower_left_center[0], self.y3d_lower_left_center[0], self.z3d_lower_left_center[0])
            elif side == 'lower-right':
                return (self.x3d_lower_right_center[0], self.y3d_lower_right_center[0], self.z3d_lower_right_center[0])
        raise ValueError('Invalid side or feature.', side, feature)


class WingSectionGeometry:
    '''
    Geometry description of the wing section (rib).
    
    Attributes
    ---------------
    parameters: dict
        dictionary containing the geometry parameters.
    fu, fl: callable
        functions to get the y coordinates of the upper and lower surfaces of the airfoil,
        given x coordinates.
    
    '''
    def __init__(self):
        
        self.parameters = {
            'xLE': 0.0, 'yLE': 0.0, 'zLE': 0.0, # leading edge coordinates (mm)
            'chord': 1.0, # chord length (mm)
            'twist': 0.0, # twist angle (degrees)
            'tmax_airfoil': 0.12, # specified airfoil maximum thickness (relative to chord)
            'xx': np.array([]), # x coordinates of the unit airfoil
            'yu': np.array([]), # y coordinates of the upper surface of the unit airfoil
            'yl': np.array([]), # y coordinates of the lower surface of the unit airfoil
            
            't_cover': 0.02, # wing cover thickness (mm)
            
            'spars': [], # spar geometries
            'stringers': [], # stringer geometries
            'cutouts': [], # cutout geometries
        }
        
        self._n_points = 101
        self.spars = []
        self.stringers = []
        self.cutouts = []
    
    def set_parameters(self, parameters):
        '''
        Set the geometry parameters of the wing section.
        
        Parameters
        ---------------
        parameters: dict
            dictionary containing the geometry parameters. 
            The keys should be the same as those in self.parameters.
        '''
        for key in parameters:
            
            if key == 'airfoil':
                # read airfoil coordinates from file
                airfoil = read_airfoil(parameters['airfoil'])
                upper, lower = split_upper_lower(airfoil)
                xx, yu, yl = reconstruct_airfoil(upper, lower, n_points=self._n_points)
                self.parameters['xx'] = xx
                self.parameters['yu'] = yu
                self.parameters['yl'] = yl
            
            if key in self.parameters:
                self.parameters[key] = parameters[key]
                
        for params in self.parameters['spars']:
            self.spars.append(SparGeometry(x=params['x'], t=params['t']))
        
        for params in self.parameters['stringers']:
            self.stringers.append(StringerGeometry(x=params['x'], t=params['t'], h=params['h'], w=params['w']))
        
        for params in self.parameters['cutouts']:
            self.cutouts.append(CutoutGeometry(x=params['x'], y=params['y'],
                                    lx=params['lx'], ly=params['ly'], r_fillet=params['r_fillet']))

        self._check_parameters()
        self._update_attributes()
        self._check_intersections()
    
    def _check_parameters(self):
        '''
        Check the validity of the geometry parameters.
        '''
        if self.xx.size == 0 or self.yu.size == 0 or self.yl.size == 0:
            raise ValueError('Airfoil coordinates are not set.')
        
        if self.xx.shape != self.yu.shape or self.xx.shape != self.yl.shape:
            raise ValueError('Inconsistent airfoil coordinates: xx, yu and yl should have the same shape.')
        
        if np.any(self.xx < 0) or np.any(self.xx > 1):
            raise ValueError('Invalid airfoil coordinates: x should be in [0,1].')
        
        if len(self.spars) < 2:
            raise ValueError('At least 2 spars are required for a valid wing box geometry.')
        else:
            for spar in self.spars:
                if not isinstance(spar, SparGeometry):
                    raise TypeError('Invalid spar geometry: each spar should be an instance of SparGeometry.')
                if spar.x <= 0 or spar.x >= 1:
                    raise ValueError('Invalid spar geometry: x should be in [0,1].')
            # sort spars by x coordinate
            self.parameters['spars'] = sorted(self.spars, key=lambda spar: spar.x)
        
        for stringer in self.stringers:
            if not isinstance(stringer, StringerGeometry):
                raise TypeError('Invalid stringer geometry: each stringer should be an instance of StringerGeometry.')
            if stringer.x <= self.x_front_spar or stringer.x >= self.x_rear_spar:
                raise ValueError('Invalid stringer geometry: x should be in [%f, %f].' % (self.x_front_spar, self.x_rear_spar))
        # sort stringers by x coordinate
        if len(self.stringers) > 0:
            self.parameters['stringers'] = sorted(self.stringers, key=lambda stringer: stringer.x)
            
        for cutout in self.cutouts:
            if not isinstance(cutout, CutoutGeometry):
                raise TypeError('Invalid cutout geometry: each cutout should be an instance of CutoutGeometry.')
            if cutout.x - cutout.lx*0.5 <= self.x_front_spar or cutout.x + cutout.lx*0.5 >= self.x_rear_spar:
                raise ValueError('Invalid cutout geometry (x): cutout should be within the front and rear spars in the x direction.')
        # sort cutouts by x coordinate
        if len(self.cutouts) > 0:
            self.parameters['cutouts'] = sorted(self.cutouts, key=lambda cutout: cutout.x)

    def _update_attributes(self):
        '''
        Update the attributes of the wing section geometry based on the parameters.
        '''
        # Convert airfoil coordinates to numpy arrays if they are not already
        self.parameters['xx'] = np.array(self.parameters['xx'])
        self.parameters['yu'] = np.array(self.parameters['yu'])
        self.parameters['yl'] = np.array(self.parameters['yl'])
        
        # Create interpolation functions for the upper and lower surfaces of the airfoil
        self.fu = interp1d(self.xx, self.yu, kind='cubic')
        self.fl = interp1d(self.xx, self.yl, kind='cubic')
        
        # Create reference x coordinates distributed between the front and rear spars
        self._xx_reference = np.linspace(0, 1, num=self._n_points, endpoint=True)
        
        # Update the attributes of the spars, stringers and cutouts
        for spar in self.spars:
            y_upper = self.fu(spar.x)
            y_lower = self.fl(spar.x)
            spar.update_attributes(y_upper, y_lower,
                                chord=self.chord, twist_degrees=self.twist_degrees,
                                xLE=self.xLE, yLE=self.yLE, zLE=self.zLE)
        for stringer in self.stringers:
            y_upper = self.fu(stringer.x)
            y_lower = self.fl(stringer.x)
            stringer.update_attributes(y_upper, y_lower,
                                chord=self.chord, twist_degrees=self.twist_degrees,
                                xLE=self.xLE, yLE=self.yLE, zLE=self.zLE)
        for cutout in self.cutouts:
            cutout.update_attributes(y_upper=None, y_lower=None,
                                chord=self.chord, twist_degrees=self.twist_degrees,
                                xLE=self.xLE, yLE=self.yLE, zLE=self.zLE)
            
        # Calculate the wing cover coordinates
        _xx_wing_cover = np.linspace(self.x_front_spar, self.x_rear_spar, num=self._n_points, endpoint=True)
        self.x3d_upper_cover, self.y3d_upper_cover = self.get_airfoil_upper_surface(
                _xx_wing_cover, input_relative_coordinates=True)
        self.x3d_lower_cover, self.y3d_lower_cover = self.get_airfoil_lower_surface(
                _xx_wing_cover, input_relative_coordinates=True)
        self.z3d_upper_cover = np.ones_like(self.x3d_upper_cover) * self.zLE
        self.z3d_lower_cover = np.ones_like(self.x3d_lower_cover) * self.zLE

    def _check_intersections(self):
        '''
        Check the intersections between the spars, covers and cutouts.
        '''
        for cutout in self.cutouts:
            hlx = cutout.lx*0.5
            hly = cutout.ly*0.5
            yu0 = self.fu(cutout.x-hlx)
            yl0 = self.fl(cutout.x-hlx)
            yu1 = self.fu(cutout.x+hlx)
            yl1 = self.fl(cutout.x+hlx)
            if cutout.y + hly >= yu0 or cutout.y - hly <= yl0 or cutout.y + hly >= yu1 or cutout.y - hly <= yl1:
                raise ValueError('Invalid cutout geometry (y): cutout should not intersect with the airfoil surfaces.')

    @property
    def xLE(self):
        '''
        x coordinate of the leading edge (mm)
        '''
        return float(self.parameters['xLE'])
    
    @property
    def yLE(self):
        '''
        y coordinate of the leading edge (mm)
        '''
        return float(self.parameters['yLE'])
    
    @property
    def zLE(self):
        '''
        z coordinate of the leading edge (mm)
        '''
        return float(self.parameters['zLE'])
    
    @property
    def chord(self):
        '''
        chord length (mm)
        '''
        return float(self.parameters['chord'])
    
    @property
    def twist_degrees(self):
        '''
        twist angle (degrees)
        '''
        return float(self.parameters['twist'])
    
    @property
    def twist_radians(self):
        '''
        twist angle (radians)
        '''
        return np.radians(self.parameters['twist'])
    
    @property
    def tmax_airfoil(self):
        '''
        specified airfoil maximum thickness (relative to chord)
        '''
        return float(self.parameters['tmax_airfoil'])

    @property
    def xx(self):
        '''
        x coordinates of the unit airfoil (np.ndarray)
        '''
        return np.atleast_1d(self.parameters['xx'])
    
    @property
    def yu(self):
        '''
        y coordinates of the upper surface of the unit airfoil (np.ndarray)
        '''
        return np.atleast_1d(self.parameters['yu'])

    @property
    def yl(self):   
        '''
        y coordinates of the lower surface of the unit airfoil (np.ndarray)
        '''
        return np.atleast_1d(self.parameters['yl'])
    
    @property
    def xx_reference(self):
        '''
        reference x coordinates (relative to chord) distributed between the front and rear spars.
        '''
        return self._xx_reference
    
    @property
    def t_cover(self):
        '''
        wing cover thickness (mm)
        '''
        return float(self.parameters['t_cover'])
    
    @property
    def n_spars(self):
        '''
        number of spars
        '''
        return len(self.spars)
    
    @property
    def x_front_spar(self):
        '''
        x coordinate of the front spar (relative to chord)
        '''
        return float(self.spars[0].x)
    
    @property
    def x_rear_spar(self):
        '''
        x coordinate of the rear spar (relative to chord)
        '''
        return float(self.spars[-1].x)

    @property
    def n_stringers(self):
        '''
        number of (pairs of) stringers
        '''
        return len(self.stringers)

    @property
    def n_cutouts(self):
        '''
        number of cutouts
        '''
        return len(self.cutouts)
    
    
    def get_airfoil_upper_surface(self, x=None, input_relative_coordinates=True):
        '''
        Get the y coordinates of the upper surface of the airfoil, given x coordinates.
        
        Parameters
        ---------------
        x: array-like or None
            x coordinates (relative to chord if input_relative_coordinates is True, otherwise in mm)
            If None, use self.xx as default x coordinates.
        
        input_relative_coordinates: bool
            whether the input x coordinates are relative to chord (True) or in mm (False).
        
        Returns
        ---------------
        new_x, new_y: ndarray
            x, y coordinates of the upper surface corresponding to the input x coordinates.
        '''
        if x is None:
            x = self.xx
        else:
            x = np.atleast_1d(x)
        
        if not input_relative_coordinates:
            x = np.array(x) / self.chord
        
        yu = self.fu(x)
        new_x, new_y = transform_curve(x, yu, scale=self.chord, rotation=self.twist_degrees, dx=self.xLE, dy=self.yLE, center=[0,0])
    
        return new_x, new_y
    
    def get_airfoil_lower_surface(self, x=None, input_relative_coordinates=True):
        '''
        Get the y coordinates of the lower surface of the airfoil, given x coordinates.
        
        Parameters
        ---------------
        x: array-like or None
            x coordinates (relative to chord if input_relative_coordinates is True, otherwise in mm)
            If None, use self.xx as default x coordinates.
        
        input_relative_coordinates: bool
            whether the input x coordinates are relative to chord (True) or in mm (False).
        
        Returns
        ---------------
        new_x, new_y: ndarray
            x, y coordinates of the lower surface corresponding to the input x coordinates.
        '''
        if x is None:
            x = self.xx
        else:
            x = np.atleast_1d(x)
        
        if not input_relative_coordinates:
            x = np.array(x) / self.chord
        
        yl = self.fl(x)
        new_x, new_y = transform_curve(x, yl, scale=self.chord, rotation=self.twist_degrees, dx=self.xLE, dy=self.yLE, center=[0,0])
    
        return new_x, new_y

    
    def get_selection_points(self, feature='cover', side='upper', index=0):
        '''
        Get the selection points of all features of the wing section.
        
        Parameters
        ---------------
        feature: str
            'cover', 'spar', 'stringer' or 'cutout',
            indicating the feature type.
        side: str
            'upper' or 'lower' for cover and stringer;
            None for spar and cutout.
        index: int
            index of the feature instance (starting from 0),
            in the x increasing order for spars, stringers, and cutouts.
        
        Returns
        ---------------
        points: list of tuples
            List of (x, y, z) coordinates (mm) of the selection points of the features.
        '''
        points = []
        if feature == 'spar':
            points.append(self.spars[index].get_selection_point(feature='spar', side=None))
        
        elif feature == 'stringer':
            points.append(self.stringers[index].get_selection_point(feature='web', side=side))
            points.append(self.stringers[index].get_selection_point(feature='flange', side=side))
        
        elif feature == 'cutout':
            for _side in ['upper', 'lower', 'left', 'right']:
                points.append(self.cutouts[index].get_selection_point(feature='edge', side=_side))
            for _side in ['upper-left', 'upper-right', 'lower-left', 'lower-right']:
                points.append(self.cutouts[index].get_selection_point(feature='fillet-curve', side=_side))

        elif feature == 'cover':
            # The cover is splitted into pieces by spars and stringers,
            # so we take the mid points of the cover segments between the spars and stringers as selection points.
            x_segments = [spar.x for spar in self.spars] + [stringer.x for stringer in self.stringers]
            x_segments.sort()
            x_mid_points = []
            for i in range(len(x_segments)-1):
                x_mid_points.append(0.5 * (x_segments[i] + x_segments[i+1]))

            if side == 'upper':
                new_x, new_y = self.get_airfoil_upper_surface(x=x_mid_points, input_relative_coordinates=True)
            elif side == 'lower':
                new_x, new_y = self.get_airfoil_lower_surface(x=x_mid_points, input_relative_coordinates=True)
            else:
                raise ValueError('Invalid side specified for cover: should be "upper" or "lower".', side)
        
            new_z = np.ones_like(new_x) * self.zLE
            points = np.concatenate((new_x[:, np.newaxis], new_y[:, np.newaxis], new_z[:, np.newaxis]), axis=1)
            points = points.tolist()
            points = [(float(p[0]), float(p[1]), float(p[2])) for p in points]
        
        else:
            raise ValueError('Invalid feature or side specified.', feature, side, index)

        return points


def plot_wing_section_geometry(wing_section_geometry):
    '''
    Plot the wing section geometry using matplotlib.
    
    Parameters
    ---------------
    wing_section_geometry: WingSectionGeometry
        the wing section geometry to be plotted.
    '''
    import matplotlib.pyplot as plt

    wsg = wing_section_geometry
    fig, ax = plt.subplots(figsize=(12, 5))

    # Airfoil outline
    xu, yu = wsg.get_airfoil_upper_surface()
    xl, yl = wsg.get_airfoil_lower_surface()
    ax.plot(xu, yu, 'k-', lw=wsg.t_cover, label='Airfoil')
    ax.plot(xl, yl, 'k-', lw=wsg.t_cover)
    ax.plot([xu[-1], xl[-1]], [yu[-1], yl[-1]], 'k-', lw=wsg.t_cover)  # trailing edge

    # Spars
    for i, spar in enumerate(wsg.spars):
        ax.plot(spar.x3d, spar.y3d, 'b-', lw=spar.t,
                label='Spar' if i == 0 else '_nolegend_')

    # Stringers (upper and lower, web and flange)
    for i, st in enumerate(wsg.stringers):
        lbl = 'Stringer' if i == 0 else '_nolegend_'
        ax.plot(st.x3d_upper_web,    st.y3d_upper_web,    'r-', lw=st.t, label=lbl)
        ax.plot(st.x3d_upper_flange, st.y3d_upper_flange, 'r-', lw=st.t, label='_nolegend_')
        ax.plot(st.x3d_lower_web,    st.y3d_lower_web,    'r-', lw=st.t, label='_nolegend_')
        ax.plot(st.x3d_lower_flange, st.y3d_lower_flange, 'r-', lw=st.t, label='_nolegend_')

    # Cutouts - assemble closed outline: left -> upper-left fillet -> upper -> upper-right fillet
    #           -> right (reversed) -> lower-right fillet -> lower (reversed) -> lower-left fillet
    for i, co in enumerate(wsg.cutouts):
        lbl = 'Cutout' if i == 0 else '_nolegend_'
        xx = np.concatenate([
            co.x3d_left_edge,
            co.x3d_upper_left_fillet[::-1],
            co.x3d_upper_edge,
            co.x3d_upper_right_fillet[::-1],
            co.x3d_right_edge[::-1],
            co.x3d_lower_right_fillet[::-1],
            co.x3d_lower_edge[::-1],
            co.x3d_lower_left_fillet[::-1],
            co.x3d_left_edge[:1],
        ])
        yy = np.concatenate([
            co.y3d_left_edge,
            co.y3d_upper_left_fillet[::-1],
            co.y3d_upper_edge,
            co.y3d_upper_right_fillet[::-1],
            co.y3d_right_edge[::-1],
            co.y3d_lower_right_fillet[::-1],
            co.y3d_lower_edge[::-1],
            co.y3d_lower_left_fillet[::-1],
            co.y3d_left_edge[:1],
        ])
        ax.plot(xx, yy, 'g-', lw=1.5, label=lbl)

    ax.set_aspect('equal')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title('Wing Section Geometry  (chord={:.0f} mm, z={:.0f} mm)'.format(
        wsg.chord, wsg.zLE))
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    return fig, ax

def plot_selection_points(wing_section_geometry, fig, ax):
    '''
    Overlay the selection points of all structural features on an existing
    wing-section plot.  Intended to be called after plot_wing_section_geometry.

    Parameters
    ---------------
    wing_section_geometry: WingSectionGeometry
        the wing section geometry whose selection points are to be plotted.
    fig, ax: matplotlib Figure and Axes
        the axes returned by plot_wing_section_geometry.

    Returns
    ---------------
    fig, ax: matplotlib Figure and Axes
    '''
    wsg = wing_section_geometry

    # --- Spars ---
    spar_configs = [
        ('spar',  'upper', 'o', 'Spar mid'),
        ('root',  'upper', '^', 'Spar root (upper)'),
        ('root',  'lower', 'v', 'Spar root (lower)'),
    ]
    for feature, side, marker, label in spar_configs:
        pts = [spar.get_selection_point(feature, side) for spar in wsg.spars]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts],
                   s=50, marker=marker, c='royalblue', zorder=6, label=label)

    # --- Stringers ---
    stringer_configs = [
        ('root',   'upper', '^', 'Stringer root'),
        ('web',    'upper', 'o', 'Stringer web mid'),
        ('corner', 'upper', 's', 'Stringer corner'),
        ('flange', 'upper', 'D', 'Stringer flange mid'),
        ('tip',    'upper', '*', 'Stringer tip'),
        ('root',   'lower', 'v', None),
        ('web',    'lower', 'o', None),
        ('corner', 'lower', 's', None),
        ('flange', 'lower', 'D', None),
        ('tip',    'lower', '*', None),
    ]
    for feature, side, marker, label in stringer_configs:
        pts = [st.get_selection_point(feature, side) for st in wsg.stringers]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts],
                   s=40, marker=marker, c='tomato', zorder=6,
                   label=label if label is not None else '_nolegend_')

    # --- Cutouts ---
    edge_sides   = ['upper', 'lower', 'left', 'right']
    corner_sides = ['upper-left', 'upper-right', 'lower-left', 'lower-right',
                    'left-upper', 'left-lower', 'right-upper', 'right-lower']
    fillet_sides = ['upper-left', 'upper-right', 'lower-left', 'lower-right']

    for j, co in enumerate(wsg.cutouts):
        for k, side in enumerate(edge_sides):
            p = co.get_selection_point('edge', side)
            ax.scatter(p[0], p[1], s=40, marker='o', c='seagreen', zorder=6,
                       label='Cutout edge' if (j == 0 and k == 0) else '_nolegend_')

        for k, side in enumerate(corner_sides):
            p = co.get_selection_point('corner', side)
            ax.scatter(p[0], p[1], s=40, marker='s', c='darkorange', zorder=6,
                       label='Cutout corner' if (j == 0 and k == 0) else '_nolegend_')

        for k, side in enumerate(fillet_sides):
            p = co.get_selection_point('fillet-curve', side)
            ax.scatter(p[0], p[1], s=40, marker='^', c='seagreen', zorder=6,
                       label='Cutout fillet curve' if (j == 0 and k == 0) else '_nolegend_')

        for k, side in enumerate(fillet_sides):
            p = co.get_selection_point('fillet-center', side)
            ax.scatter(p[0], p[1], s=30, marker='x', c='seagreen', zorder=6,
                       label='Cutout fillet center' if (j == 0 and k == 0) else '_nolegend_')

    # --- Cover ---
    for k, side in enumerate(['upper', 'lower']):
        pts = wsg.get_selection_points('cover', side=side)
        ax.scatter([p[0] for p in pts], [p[1] for p in pts],
                   s=50, marker='P', c='mediumpurple', zorder=6,
                   label='Cover (%s)'%(side) if k == 0 else 'Cover (%s)'%(side))

    ax.legend(loc='upper right', fontsize=7, ncol=2)
    return fig, ax


if __name__ == '__main__':
    import json
    import matplotlib.pyplot as plt

    path = os.path.dirname(os.path.abspath(__file__))
    fname = os.path.join(path, 'default-parameters.json')
    with open(fname, 'r') as f:
        parameters = json.load(f)

    pGeo = parameters['pGeo']
    pMesh = parameters['pMesh']
    pRun = parameters['pRun']

    for section_params in pGeo['sections']:
        section_params['airfoil'] = os.path.join(path, section_params['airfoil'])
        wsg = WingSectionGeometry()
        wsg.set_parameters(section_params)
        fig, ax = plot_wing_section_geometry(wsg)
        plot_selection_points(wsg, fig, ax)
        plt.tight_layout()
        plt.show()
    