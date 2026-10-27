'''
Utility functions for defining composite layups in Abaqus.
'''

from AbaqusTools import IS_ABAQUS

if IS_ABAQUS:

    from abaqus import *
    from abaqusConstants import *
    from symbolicConstants import *


def create_shell_CompositeLayup_of_set(myPrt, name_set, ply_thickness, ply_angles,
                name_surface, primaryAxisVector,
                symmetric=True, numIntPoints=3, material_name='IM7/8551-7'):
    '''
    Create composite layup of shell elements for a set of faces.
    - Discrete layup orientation
    '''
    if name_set not in myPrt.sets.keys():
        raise Exception('Set %s not found in part %s!'%(name_set, myPrt.name))
    
    if name_surface not in myPrt.surfaces.keys():
        raise Exception('Surface %s not found in part %s!'%(name_surface, myPrt.name))
    
    compositeLayup = myPrt.CompositeLayup(name=name_set, description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=symmetric, 
        thicknessAssignment=FROM_SECTION)
    
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
                useDensity=OFF)
    
    compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
        additionalRotationType=ROTATION_NONE, angle=0.0, 
        additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
        normalAxisDefinition=SURFACE,
        normalAxisRegion= myPrt.surfaces[name_surface], 
        normalAxisDirection=AXIS_3, flipNormalDirection=False, 
        primaryAxisDefinition=VECTOR,
        primaryAxisVector=primaryAxisVector,
        primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)

    #* Plies
    for i in range(len(ply_angles)):
    
        compositeLayup.CompositePly(suppressed=False, plyName='ply-%d'%(i+1), 
            region=myPrt.sets[name_set], material=material_name, 
            thicknessType=SPECIFY_THICKNESS, thickness=ply_thickness, 
            orientationType=SPECIFY_ORIENT, orientationValue=ply_angles[i], 
            additionalRotationType=ROTATION_NONE,   # CSYS in `Plies` for rotation angle
            additionalRotationField='',             # CSYS in `Plies` for rotation angle
            axis=AXIS_3, angle=0.0,                 # CSYS in `Plies` for rotation angle
            numIntPoints=numIntPoints)

