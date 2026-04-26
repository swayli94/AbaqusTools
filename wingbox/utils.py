'''

'''
import numpy as np
from scipy.interpolate import interp1d


def read_airfoil(filename):
    '''
    Read airfoil coordinates from file.
    
    Parameters
    -------------
    filename: str
        file name
    
    Returns
    -------------
    airfoil: ndarray [:,2]
        airfoil coordinates. It has a unit chord length, x ranges in [0,1].
    '''
    airfoil = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            line = line.split()
            airfoil.append([float(line[0]), float(line[1])])
    
    return np.array(airfoil)

def split_upper_lower(airfoil, remove_trailing_edge=False):
    '''
    Split the airfoil into upper and lower surfaces.
    
    Parameters
    ---------------
    airfoil: ndarray [:,2]
        airfoil coordinates
        
    remove_trailing_edge: bool
        whether remove trailing edge thickness.
        If so, the airfoil thickness will change.
    
    Returns
    ----------------
    upper, lower: ndarray [:,2]
        coordinates of the upper and lower airfoil surfaces.
        Both contains the leading edge point.
    '''
    index_LE = np.argmin(airfoil[:,0])
    upper = airfoil[:index_LE+1,:]; upper = upper[::-1,:]
    lower = airfoil[index_LE:,:]
    
    if remove_trailing_edge:
        for i in range(upper.shape[0]):
            upper[i,1] -= upper[i,0]*upper[-1,1]
        for i in range(lower.shape[0]):
            lower[i,1] -= lower[i,0]*lower[-1,1]
            
    return upper, lower

def dist_clustcos(nn, a0=0.0079, a1=0.96, beta=1.0):
    '''
    Point distribution on x-axis [0, 1]. (More points at both ends)

    Parameters
    ----------
    nn: int
        total amount of points
        
    a0: float
        Parameter for distributing points near x=0.
        Smaller a0, more points near x=0.
        
    a1: float
        Parameter for distributing points near x=1.
        Larger a1, more points near x=1.
        
    beta: float
        Parameter for distribution points.
    
    Examples
    ---------
    >>> xx = dist_clustcos(n, a0, a1, beta)

    '''
    aa = np.power((1-np.cos(a0*np.pi))/2.0, beta)
    dd = np.power((1-np.cos(a1*np.pi))/2.0, beta) - aa
    yt = np.linspace(0.0, 1.0, num=nn)
    a  = np.pi*(a0*(1-yt)+a1*yt)
    xx = (np.power((1-np.cos(a))/2.0,beta)-aa)/dd

    return xx

def reconstruct_airfoil(upper, lower, n_points=101):
    '''
    Reconstruct the airfoil coordinates from upper and lower surfaces.
    
    Parameters
    ---------------
    upper, lower: ndarray [:,2]
        coordinates of the upper and lower airfoil surfaces.
        Both contains the leading edge point.
    
    Returns
    ----------------
    xx, yu, yl: ndarray
        x coordinates, y coordinates of upper and lower surfaces.
    '''
    xx = dist_clustcos(n_points, a0=0.0079, a1=0.96, beta=1.0)
    f_upper = interp1d(upper[:,0], upper[:,1], kind='cubic')
    f_lower = interp1d(lower[:,0], lower[:,1], kind='cubic')
    yu = f_upper(xx)
    yl = f_lower(xx)
    return xx, yu, yl


def norm_vec(vector):
    return np.array(vector)/np.linalg.norm(vector)

def triangle_area(x1, y1, x2, y2, x3, y3):
    '''
    Calculate the area of a triangle specified by three points
    '''
    return 0.5*abs((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))

def rotate3d(points, alpha, beta, gamma, origin=[0, 0, 0]):
    '''
    Rotate points about a given point (origin).
    
    Parameters
    ------------------
    points: ndarray [:,3]
        points to be rotated
    
    alpha: float
        first rotate alpha about the x-axis
    
    beta: float
        then rotate beta about the y-axis
    
    gamma: float
        at last, rotate gamma about the z-axis
        
    origin: list [3]
        the origin of rotation
        
    Returns
    ----------------
    new_points: ndarray [:,3]
        rotated points
    '''
    # https://zhuanlan.zhihu.com/p/388164543
    
    alpha, beta, gamma = np.deg2rad(alpha), np.deg2rad(beta), np.deg2rad(gamma)
    
    Rx = [[            1,              0,             0,                                                       0],
          [            0,  np.cos(alpha),-np.sin(alpha), origin[1]*(1 - np.cos(alpha)) + origin[2]*np.sin(alpha)],
          [            0,  np.sin(alpha), np.cos(alpha), origin[2]*(1 - np.cos(alpha)) - origin[1]*np.sin(alpha)],
          [            0,              0,             0,                                                       1]]
    Ry = [[ np.cos(beta),              0,  np.sin(beta), origin[0]*(1 - np.cos(beta))  - origin[2]*np.sin(beta) ],
          [            0,              1,             0,                                                       0],
          [-np.sin(beta),              0,  np.cos(beta), origin[2]*(1 - np.cos(beta))  + origin[0]*np.sin(beta) ],
          [            0,              0,             0,                                                       1]]
    Rz = [[np.cos(gamma), -np.sin(gamma),             0, origin[0]*(1 - np.cos(gamma)) + origin[1]*np.sin(gamma)],
          [np.sin(gamma),  np.cos(gamma),             0, origin[1]*(1 - np.cos(gamma)) - origin[0]*np.sin(gamma)],
          [            0,              0,             1,                                                       0],
          [            0,              0,             0,                                                       1]]
    
    Rx = np.array(Rx); Ry = np.array(Ry); Rz = np.array(Rz)
    
    Pp = np.array([points[:,0], points[:,1], points[:,2], [1.0]*points.shape[0]])
    
    return np.dot(np.dot(np.dot(Rx, Ry), Rz), Pp)


def transform_curve(x, y, scale=1.0, rotation=0.0, dx=0.0, dy=0.0, center=[0,0]):
    '''
    Transform a curve, in the order of scaling, rotation and translation.
    
    Parameters
    -------------------
    x, y: ndarray
        x and y coordinates of the curve
        
    scale: float
        scaling factor
        
    dx, dy: float
        translation distance in x, y directions
        
    rotation: float
        rotation angle (degree)
        
    center: array [2]
        the center point for scaling and rotation
    
    Returns
    ------------
    new_x, new_y: ndarray
        new curve coordinates
    '''
    new_x = (np.array(x) - center[0])*scale
    new_y = (np.array(y) - center[1])*scale
    
    n = len(x)
    if n != len(y):
        raise ValueError('x and y should have the same length.')

    angle = np.radians(rotation)
    cc = np.cos(angle)
    ss = np.sin(angle)
    for i in range(n):
        new_x[i] = new_x[i]*cc - new_y[i]*ss
        new_y[i] = new_x[i]*ss + new_y[i]*cc

    new_x += center[0] + dx
    new_y += center[1] + dy

    return new_x, new_y

