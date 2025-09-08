'''
This is a python implementation of LaRC05 failure criterion.
It is translated from the user subroutine Fortran codes for Abaqus, including:

    1. uvarm.f90
    2. failureCriteria.f90
    3. plyTools.f90 

Reference:

    LaRC05 developed by Miguel A.S. Matos and Silvestre T. Pinho @ 17.04.2019
    Imperial College London, Aeronautics. 
'''

import numpy as np

from typing import Final, Any, Tuple


class UVARM(object):
    '''
    User-defined output variables.
    
    https://docs.software.vt.edu/abaqusv2023/English/?show=SIMACAESUBRefMap/simasub-c-uvarm.htm
    
    
    Parameters
    ------------------        
    nSCply: int

        Number of stress components at this point, 3 or 6.

    material: str
    
        Material name, e.g., 'IM7/8551-7'.


    Attributes
    ------------------
    NFI: int
    
        Number of failure indexes, default is 5. Including:
        
        1) Matrix cracking;
        2) Matrix splitting;
        3) Fibre tension;
        4) Fibre kinking;
        5) Transverse inter-bundle failure mode;
    
    NUVARM: int
    
        Number of user-defined output variables, default is 7. Including:

        1-5) Failure indexes;
        6) Maximum failure index;
        7) Failure mode;
    
    LIMIT_UVARM: float

        Constant. Maximum value of the max failure index (UVARM6), default is 1.0.

    isElement3D: bool
    
        Whether processing 3D elements or not.

    '''
    def __init__(self, nSCply: int, material='IM7/8551-7') -> None:
        
        #* Constant parameters
        self.NFI         : Final[int] = 5
        self.LIMIT_UVARM : Final[float] = 1.0
        
        #* Global parameters
        self.nSCply = nSCply
        self.NUVARM = self.NFI + 2
        
        if self.nSCply == 3:
            self.isElement3D = False
        elif self.nSCply == 6:
            self.isElement3D = True
        else:
            print()
            print('Error [UVARM]: __init__')
            print('    Wrong number of stress components input [nSCply]: ', nSCply)
            print()
            raise Exception

        #* Initialize ply properties
        rProps = self.get_property(material)
        
        self.ply = PlyProperty(rProps)
        
        self.FISolver = FailureCriteria(self.nSCply)

    def get_uvarm(self, plyStresses: np.ndarray, oldUVARM=None, limitFIDen=False) -> np.ndarray:
        '''
        Calculate the user-defined output variable based on the in-situ stresses and material properties.
        
        Parameters
        ------------------
        plyStresses: ndarray [nSCply]
        
            Ply stresses in local reference frame
        
        oldUVARM: ndarray [7], or None
        
            History UVARM of this ply. If None, then zeros.

        limitFIDen: bool
        
            Flag for denominator limitation for damage propagation

        Returns
        ------------------
        currentUVARM: ndarray [7]
        
            An array containing the user-defined output variables. 
            In Abaqus, these are passed in as the values at the beginning of the increment and 
            must be returned as the values at the end of the increment.
        
        '''        
        if oldUVARM is None:
            oldUVARM = np.zeros(self.NUVARM)
        
        plyIndexes = oldUVARM[:self.NFI].copy()
        
        #* Calculate failure indexes
        plyIndexes = self.FISolver.completeCriteria(self.ply, plyStresses, plyIndexes, limitFIDen)
        
        #* Assemble UVARM
        currentUVARM = np.zeros_like(oldUVARM)
        currentUVARM[:self.NFI] = plyIndexes
        
        # UVARM 6 -> maximum index
        currentUVARM[5] = min( np.max(plyIndexes), self.LIMIT_UVARM )
        
        # UVARM 7 -> failure mode
        if currentUVARM[5] >= self.LIMIT_UVARM - 1e-2:
            
            if oldUVARM[6] > 0.0:
                currentUVARM[6] = oldUVARM[6]
            else:
                currentUVARM[6] = np.argmax(plyIndexes) + 1.0

        return currentUVARM

    def get_property(self, material='IM7/8551-7') -> dict:
        '''
        Perform interpolation of properties defined in the property table in Abaqus/Standard.
        
        https://docs.software.vt.edu/abaqusv2023/English/?show=SIMACAESUBRefMap/simasub-c-tablecollection.htm#simasub-c-tableCollection-getPropertyTables
        
        
        Returns
        ------------------
        rProps: dict
        
            Dictionary containing material properties.
            
            Use the default value if the key is missing, or the value is None.

        '''
        rProps_IM7_8551_7 = {
            'E11':      1.6500E+5,  # Young’s modulus along the fibre direction
            'E22':      8.4000E+3,  # Young’s modulus along transverse direction
            'nu12':     3.4000E-1,  # Poisson’s ratio
            'G12':      5.6000E+3,  # In plane shear modulus 1-2 plane
            'Xt':       2.5600E+3,  # Longitudinal tensile strength
            'Xc':       1.5900E+3,  # Longitudinal compressive strength
            'Yt':       7.3000E+1,  # Transverse tensile strength 
            'Yc':       1.8500E+2,  # Transverse compressive strength
            'Sl':       9.0000E+1,  # In plane shear strength
            'a0':       5.3000E+1,  # Fracture angle for pure compression (53 degree)
            'nL':       8.2000E-2,  # Longitudinal shear friction coefficient
            'nT':       None,       # Transverse shear friction coefficient
            'ILSS':     9.0000E+1,  # Inter-laminar shear strength
            'Zt':       6.3000E+1,  # Tensile strength in the 3rd direction
            'G1cMat':   2.1000E-1,  # Matrix mode I  critical energy release rate
            'G2cMat':   8.0000E-1,  # Matrix mode II critical energy release rate
            'G1cFibT':  9.2000E+1,  # Fibre tensile  critical energy release rate
            'G1cFibK':  8.0000E+1,  # Fibre kinking  critical energy release rate
            'GAlphaM':  1.2100E+0,  # Matrix mixed-mode power law exponent (1.21)
        }
        
        if material == 'IM7/8551-7':
            rProps = rProps_IM7_8551_7
            
        else:
            print()
            print('Error [UVARM]: getProperty')
            print('    Wrong material input [material]: ', material)
            print()
            raise Exception
            
        return rProps
    

class PlyProperty(object):
    '''
    Ply's material properties.
    
    Parameters
    ------------------
    rProps: dict
    
        Dictionary containing material properties.
    
    
    Attributes
    ------------------
    E11, E22: float
        elastic moduli
        
    nu12, nu21: float
        Poisson's ratio
        
    G12: float
        shear moduli
        
    Xt, Xc, Yt, Yc: float
        strengths
        
    Sl, a0, nL, nT, St: float
        additional properties
        
    YtIS, SlIS, StIS: float
        in-situ effects
        
    phiC: float
        critical fibre misalignment angle
        
    ILSS: float
        Inter-laminar shear strength
    
    Zt: float
        Tensile strength in the 3rd direction
    
    GValues: ndarray [5]
        Failure energies
        
    G1cMat,G2cMat: float
        Matrix mode I, II critical energy release rate

    G1cFibT,G1cFibK: float
        Fibre tensile, kinking critical energy release rate
        
    GAlphaM: float
        Matrix mixed-mode power law exponent (default, 1.21)
        
    E11comp: float
        Compression modulus
        
    E33: float
        Through the thick
        
    nu23, nu13: float
        Poisson's
        
    G13, G23: float
        more shear components
        
    viscosity: float
        viscous regularization
        
    PTYP: int
        ply type
        
    TENCOMP: bool
        tension compression flag

    '''
    def __init__(self, rProps: dict, 
                    DEFAULT_A0=53.0, DEFAULT_AG=1.21, 
                    DEFAULT_NU23=0.5, DEFAULT_PTYP=1) -> None:
        '''
        Initialize the ply by properties
        
        Parameters
        ------------------
        rProps: dict
        
            Dictionary containing material properties.
            
        DEFAULT_A0: float

            Default fracture angle for pure compression (53 degrees)
            
        DEFAULT_AG: float
        
            Default matrix mixed-mode power law exponent (1.21)
            
        DEFAULT_NU23: float
        
            Default Poisson's ratio in 23-plane (0.5)
            
        DEFAULT_PTYP: int

            Default ply type (1 means UD)
        '''
        PLY_TYPE_UD             : Final[int] = 1
        PLY_TYPE_THICK_EMBEDDED : Final[int] = 2
        PLY_TYPE_THIN_EMBEDDED  : Final[int] = 3
        PLY_TYPE_THIN_OUTER     : Final[int] = 4
        PLY_TYPE_EMBEDDED_AUTO  : Final[int] = 5
        
        #* ---------------------------------------------------
        #* Load properties
        #* ---------------------------------------------------
        #* Basic properties
        self.E11 =  rProps['E11']
        self.E22 =  rProps['E22']
        self.nu12 = rProps['nu12']
        self.G12 =  rProps['G12']
        self.Xt =   rProps['Xt']
        self.Xc =   rProps['Xc']
        self.Yt =   rProps['Yt']
        self.Yc =   rProps['Yc']
        self.Sl =   rProps['Sl']
        
        self.a0 =   self.get_value(rProps, 'a0', default= DEFAULT_A0) / 180.0 * np.pi # radian

        self.nL =   self.get_value(rProps, 'nL', 
                                    default= - self.Sl*np.cos(2*self.a0) / (self.Yc*np.cos(self.a0)**2))
        
        self.nT =   self.get_value(rProps, 'nT', default= -1/np.tan(2*self.a0))

        #* Failure energies
        self.G1cMat  = self.get_value(rProps, 'G1cMat',  None)
        self.G2cMat  = self.get_value(rProps, 'G2cMat',  None)
        self.G1cFibT = self.get_value(rProps, 'G1cFibT', None)
        self.G1cFibK = self.get_value(rProps, 'G1cFibK', None)
        self.GAlphaM = self.get_value(rProps, 'GAlphaM', DEFAULT_AG)
        self.GValues = np.array([self.G1cMat, self.G2cMat, self.G1cFibT, self.G1cFibK, self.GAlphaM])
        
        #* Compression modulus
        self.TENCOMP = self.get_value(rProps, 'TENCOMP', default= True)
        self.E11comp = self.get_value(rProps, 'E11comp', default= self.E11)

        #* 3D properties
        self.E33    = self.get_value(rProps, 'E33',     default= self.E22)
        self.nu13   = self.get_value(rProps, 'nu13',    default= self.nu12)
        self.nu23   = self.get_value(rProps, 'nu23',    default= DEFAULT_NU23)
        self.G13    = self.get_value(rProps, 'G13',     default= self.G12)
        self.G23    = self.get_value(rProps, 'G23',     default= self.E33/(2*(1+self.nu23)))

        #* 3D NCF failure indexes
        self.ILSS   = self.get_value(rProps, 'ILSS',    default= None)
        self.Zt     = self.get_value(rProps, 'Zt',      default= 0.5*self.Yt)
        
        if self.ILSS is None:
            self.Zt = None
            
        #* Viscous stabilization
        self.viscosity = self.get_value(rProps, 'viscosity', default= None)

        #* Ply type 
        self.PTYP   = self.get_value(rProps, 'PTYP', default= DEFAULT_PTYP)


        #* ---------------------------------------------------
        #* Calculate dependent properties
        #* ---------------------------------------------------
        #* Shear strength
        self.St = self.Yc * np.cos(self.a0) * (np.sin(self.a0) + np.cos(self.a0)/np.tan(2*self.a0))

        #* Poisson ratio
        if self.TENCOMP:
            aux1 = self.E11comp
        else:
            aux1 = self.E11
            
        self.nu21 = self.nu12 * self.E22 / aux1

        #* Critical fibre misalignment angle (radian)
        aux1 = self.Sl / self.Xc
        aux2 = 2.0 * (self.nL + aux1)
        aux1 = 1.0 - np.sqrt(1-2*aux1*aux2)
        self.phiC = np.arctan(aux1/aux2)

        #* In-situ effects
        if self.PTYP == PLY_TYPE_UD:
            
            self.StIS = self.St
            self.YtIS = self.Yt
            self.SlIS = self.Sl

        elif self.PTYP == PLY_TYPE_THICK_EMBEDDED:
            
            sq2 = np.sqrt(2)
            self.StIS = self.St * sq2
            self.YtIS = self.Yt * sq2 * 1.12
            self.SlIS = self.Sl * sq2
    
        elif self.PTYP == PLY_TYPE_THIN_EMBEDDED:
            
            G1c  = rProps['G1c']
            G2c  = rProps['G2c']
            Th   = rProps['Th']
            aa   = 2 * (1/self.E22 - self.nu21/self.E11)
            self.YtIS = np.sqrt( 8 * G1c / (np.pi * Th * aa) )
            self.SlIS = np.sqrt( 8 * G2c * self.G12 / (np.pi * Th) )
            self.StIS = self.SlIS

        elif self.PTYP == PLY_TYPE_THIN_OUTER:
            
            G1c  = rProps['G1c']
            G2c  = rProps['G2c']
            Th   = rProps['Th']
            aa   = 2 * (1/self.E22 - self.nu21/self.E11)
            self.YtIS = np.sqrt( 4 * G1c / (np.pi * Th * aa) )
            self.SlIS = np.sqrt( 4 * G2c * self.G12 / (np.pi * Th) )
            self.StIS = self.SlIS * np.sqrt(2)

        elif self.PTYP == PLY_TYPE_EMBEDDED_AUTO:
            
            G1c  = rProps['G1c']
            G2c  = rProps['G2c']
            Th   = rProps['Th']
            sq2  = np.sqrt(2)
            aa   = 2 * (1/self.E22 - self.nu21/self.E11)
            aux1 = self.Yt * sq2 * 1.12
            aux2 = np.sqrt( 8 * G1c / (np.pi * Th * aa) )
            self.YtIS = max( aux1, aux2 )
            aux1 = self.Sl * sq2
            aux2 = np.sqrt( 8 * G2c * self.G12 / (np.pi * Th) )
            self.SlIS = max( aux1, aux2 )
            aux1 = self.St * sq2
            self.StIS = max( aux1, aux2 )
            
        else:
            print()
            print('Error [PlyProperty]: __init__')
            print('    Wrong Ply type input: ', self.PTYP)
            print()
            raise Exception
                
    @staticmethod
    def get_value(source: dict, key: str, default=None) -> Any:
        '''
        Get value from a dictionary `source` by key.
        '''
        if key in source.keys():
            value = source[key]
        else:
            value = None
            
        if value is None:
            return default
        else:
            return value
    
    @property
    def a0_degree(self) -> float:
        '''
        Fracture angle for pure compression (degree), default `a0` is in radian.
        '''
        return self.a0/np.pi*180.0
        
    def plyEvaluateCriteria(self, SST, SSL, SN, limitFIDen=False) -> float:
        '''
        Evaluates the LaRC05 failure index, used by modes: 1) matrix cracking, 
        2) matrix splitting, and 4) fibre kinking.
        
        Parameters
        --------------------
        SST: float
        
            Transverse shear stress 
        
        SSL: float
        
            Longitudinal shear stress 
            
        SN: float
        
            Normal stress
        
        Returns
        --------------------
        FI: float

            Failure index
        
        Notes
        --------------------
        It was first added @ 18.06.2019, modified @ 13.02.2020.
        So that, the beneficial effect of normal stress is limited by its value 
        when pure compressive failure happens Yc * (cos(ALPHA0)^2).
        '''
        if limitFIDen:
            SNDEN = max( -self.Yc*np.cos(self.a0)**2, min(0.0, SN) )
        else:
            SNDEN = SN

        aux1 = SST / (self.StIS - self.nT * SNDEN)
        aux2 = SSL / (self.SlIS - self.nL * SNDEN)
        aux3 = max(SN/self.YtIS, 0.0)

        return np.sqrt(aux1**2 + aux2**2 + aux3**2)

    
class FailureCriteria(object):
    '''
    LaRC05 failure criteria.
    
    Five failure modes are evaluated. Each of the modes has an associated failure index FI 
    which has a unit value when failure occurs. 
    
    
    Parameters
    ------------------        
    nSCply: int

        Number of stress components at this point.
    
    
    Attributes
    ------------------
    n11, n22, n33, n12, n13, n23: int
    
        Constant, index of Sij in the stress vector

    NZNEG: float
    
        Constant, NZNEG = -1E-12
    
    NZSTRESS: float
    
        Constant, NZSTRESS = 1E-6
        
    NEG_STRESS_THS: float
    
        Constant, NEG_STRESS_THS = -1E-3
    
    KINK_ANGLE_SPACING: float

        Constant, KINK_ANGLE_SPACING = 30.0
    
    isElement3D: bool
    
        Whether processing 3D elements or not.
    
    '''
    def __init__(self, nSCply: int) -> None:
        
        self.n11 : Final[int] = 0
        self.n22 : Final[int] = 1
        self.n33 : Final[int] = 2
        self.n12 : Final[int] = 3
        self.n13 : Final[int] = 4
        self.n23 : Final[int] = 5
        
        self.n11_2d : Final[int] = 0
        self.n22_2d : Final[int] = 1
        self.n12_2d : Final[int] = 2
        
        self.NZSTRESS           : Final[float] = 1E-6
        self.NEG_STRESS_THS     : Final[float] = -1E-3  # Modified to tolerate near zero stresses @ 22.11.2019
        self.KINK_ANGLE_SPACING : Final[float] = 30.0
        
        self.nSCply = nSCply
        
        if self.nSCply == 3:
            self.isElement3D = False
        elif self.nSCply == 6:
            self.isElement3D = True
        else:
            print()
            print('Error [UVARM]: __init__')
            print('    Wrong number of stress components input [nSCply]: ', nSCply)
            print()
            raise Exception
        
    def rotateStress(self, inputStress: np.ndarray, angleRad: float, axis: int) -> np.ndarray:
        '''
        Rotates the stress vector around the specified axis.
        
        Parameters
        ------------------
        inputStress: ndarray [6]

            Stress vector
            
        angleRad: float
        
            Rotation angle in radian
            
        axis: int
        
            Axis (0,1,2) represents (x,y,z), respectively
        
        Returns
        ------------------
        rotS: ndarray [6]

            Rotated stress vector
        '''
        rotS = np.zeros(6)
        
        c1 = np.cos(angleRad)
        c2 = c1*c1
        s1 = np.sin(angleRad)
        s2 = s1*s1
        cs = c1*s1
        
        if axis == 0:
            
            rotS[self.n11] =     inputStress[self.n11]
            rotS[self.n22] =  c2*inputStress[self.n22] + s2*inputStress[self.n33] + 2*cs*inputStress[self.n23]
            rotS[self.n33] =  s2*inputStress[self.n22] + c2*inputStress[self.n33] - 2*cs*inputStress[self.n23]
            rotS[self.n12] =  s1*inputStress[self.n13] + c1*inputStress[self.n12]
            rotS[self.n13] =  cs*inputStress[self.n13] - s1*inputStress[self.n12]
            rotS[self.n23] = -cs*inputStress[self.n22] + cs*inputStress[self.n33] + (c2-s2)*inputStress[self.n23]
            
        elif axis == 1:
            
            rotS[self.n11] =  s2*inputStress[self.n33] + c2*inputStress[self.n11] - 2*cs*inputStress[self.n13]
            rotS[self.n22] =     inputStress[self.n22]
            rotS[self.n33] =  c2*inputStress[self.n33] + s2*inputStress[self.n11] + 2*cs*inputStress[self.n13]
            rotS[self.n12] =  cs*inputStress[self.n12] - s1*inputStress[self.n23]
            rotS[self.n13] = -cs*inputStress[self.n33] + cs*inputStress[self.n11] + (c2-s2)*inputStress[self.n13]
            rotS[self.n23] =  s1*inputStress[self.n12] + c1*inputStress[self.n23]

        elif axis == 2:
            
            rotS[self.n11] =  c2*inputStress[self.n11] + s2*inputStress[self.n22] + 2*cs*inputStress[self.n12]
            rotS[self.n22] =  s2*inputStress[self.n11] + c2*inputStress[self.n22] - 2*cs*inputStress[self.n12]
            rotS[self.n33] =     inputStress[self.n33]
            rotS[self.n12] = -cs*inputStress[self.n11] + cs*inputStress[self.n22] + (c2-s2)*inputStress[self.n12]
            rotS[self.n13] =  s1*inputStress[self.n23] + c1*inputStress[self.n13]
            rotS[self.n23] =  cs*inputStress[self.n23] - s1*inputStress[self.n13]

        else:
            print()
            print('Error [FailureCriteria]: rotateStress')
            print('    Wrong axis input: ', axis)
            print()
            raise Exception

        return rotS
    
    def get_matrix_cracking(self, ply: PlyProperty, stressVector: np.ndarray, oldFIndex: float, 
                                matrixAngles: list, limitFIDen=False) -> float:
        '''
        Get the failure index of matrix cracking.
        
        Intra-laminar matrix failure considers the transverse shear, longitudinal shear and normal stress 
        on the matrix fracture plane.
        
        Parameters
        ------------------
        ply: PlyProperty
        
            Ply material properties.
            
        stressVector: ndarray [6]
        
            3D stresses in local reference frame
            
        oldFIndex: float
        
            History failure index
            
        matrixAngles: list
        
            A list of angles for potential fracture planes
        
        limitFIDen: bool
        
            Flag for denominator limitation for damage propagation
        
        Returns
        ------------------
        FIndex: float
        
            Failure index
        '''
        currentMax = 0.0
        
        for angle in matrixAngles:
            
            c1 = np.cos(angle)
            s1 = np.sin(angle)
            c2 = np.cos(2*angle)
            s2 = np.sin(2*angle)

            # Potential fracture plane
            # The components sigma_N, tau_T, tau_L are defined in the matrix crack plane at angle 
            sn  = stressVector[self.n22] * c1**2 \
                + stressVector[self.n33] * s1**2 \
                + stressVector[self.n23] * s2
            
            sst = stressVector[self.n22] * (-0.5*s2) \
                + stressVector[self.n33] * 0.5*s2 \
                + stressVector[self.n23] * c2
            
            ssl = stressVector[self.n12] * c1 \
                + stressVector[self.n13] * s1

            # Calculate failure index
            currentMax = max(currentMax, ply.plyEvaluateCriteria(sst, ssl, sn, limitFIDen))
            
        return max(oldFIndex, currentMax)
    
    def get_fibre_tension(self, ply: PlyProperty, stressVector: np.ndarray, oldFIndex: float) -> float:
        '''
        Get the failure index of fibre tension.
        
        Occurring when the tensile stress along the fibre direction reaches XT. 
        
        Parameters
        ------------------
        ply: PlyProperty
        
            Ply material properties.
            
        stressVector: ndarray [6]
        
            3D stresses in local reference frame
            
        oldFIndex: float
        
            History failure index
                    
        Returns
        ------------------
        FIndex: float
        
            Failure index
        '''
        currentMax = max(0.0, stressVector[0]) / ply.Xt
            
        return max(oldFIndex, currentMax)
    
    def get_matrix_splitting_n_fibre_kinking(self, ply: PlyProperty, stressVector: np.ndarray, 
                        oldFI_MS: float, oldFI_FK: float, kinkAngles: list, limitFIDen=False) -> Tuple[float, float]:
        '''
        Get the failure index of matrix splitting and fibre kinking.
        
        Matrix splitting: 
        
            Shear dominated matrix failure, under longitudinal compression below 
            |sigma_11|< XC/2, sigma_11<0. 
            Failure index equation is evaluated at all possible kink/split angles psi.
            to determine the maximum index value. 
            The stresses are calculated in the fibre misalignment plane m, from the kink-band plane psi.
            The fibre misalignment angle phi is calculated from the critical phi_c angle.
        
        Fibre kinking:
        
            Assumed to result from shear-dominated matrix failure, under significant longitudinal 
            compression |sigma_11|>= XC/2, sigma_11<0.
            These stress components are respective to the fibre misalignment frame phi as described for 
            the matrix splitting failure index. 
        
        
        Parameters
        ------------------
        ply: PlyProperty
        
            Ply material properties.
            
        stressVector: ndarray [6]
        
            3D stresses in local reference frame
            
        oldFI_MS, oldFI_FK: float
        
            History failure index
            
        kinkAngles: list
        
            A list of angles for potential kink/split planes
        
        limitFIDen: bool
        
            Flag for denominator limitation for damage propagation
        
        Returns
        ------------------
        fi_MS, fi_FK: float
        
            Failure index
        '''
        fi_MS = oldFI_MS
        fi_FK = oldFI_FK
        
        if stressVector[0] >= self.NEG_STRESS_THS:
            
            return oldFI_MS, oldFI_FK
            
        currentMax = 0.0
    
        # Find psi (kink band angle)
        aux1 = (ply.G12 - ply.Xc) * ply.phiC
        
        for angle in kinkAngles:
            
            # Rotate around x, fixed 16.05.2019
            auxStresses = self.rotateStress(stressVector, angle, axis=0)
            aux2 = abs(auxStresses[3]) + aux1
            
            s1mod = max( -ply.Xc, auxStresses[0])
            s2mod = max( -ply.Yc, min(auxStresses[1], ply.Yt) )
            aux3  = ply.G12 + s1mod - s2mod

            # Fibre misalignment angle (radian)
            phi = np.sign(auxStresses[3]) * aux2/aux3
            
            # Rotate around z, fixed 16.05.2019
            auxStresses = self.rotateStress(auxStresses, phi, axis=2)
            
            currentIndex = ply.plyEvaluateCriteria(
                auxStresses[self.n23], auxStresses[self.n12], auxStresses[self.n22], limitFIDen)
        
            currentMax = max(currentMax, currentIndex)
        
        if stressVector[self.n11] < -0.5 * ply.Xc:
            # Fibre kinking
            fi_FK = max(currentMax, oldFI_FK)
        
        else:
            # Matrix splitting
            fi_MS = max(currentMax, oldFI_MS)
    
        return fi_MS, fi_FK
    
    def get_matrix_interface(self, ply: PlyProperty, stressVector: np.ndarray, oldFIndex: float) -> float:
        '''
        Get the failure index of transverse inter-bundle failure mode.
        
        Failure on a plane normal to the thickness direction, at the interface between fibre bundle and matrix.
        
        Parameters
        ------------------
        ply: PlyProperty
        
            Ply material properties.
            
        stressVector: ndarray [6]
        
            3D stresses in local reference frame
            
        oldFIndex: float
        
            History failure index
                    
        Returns
        ------------------
        FIndex: float
        
            Failure index
        '''
        if stressVector[self.n33] > self.NZSTRESS:
            
            aux1 = stressVector[self.n13] / ply.ILSS
            aux2 = stressVector[self.n23] / ply.ILSS
            aux3 = stressVector[self.n33] / ply.Zt

            currentIndex = np.sqrt(aux1**2 + aux2**2 + aux3**2)
            
        else:
            
            currentIndex = 0.0

        return max(currentIndex, oldFIndex)
        
    def completeCriteria(self, ply: PlyProperty, plyStresses: np.ndarray, oldPlyIndexes: np.ndarray, 
                            limitFIDen=False) -> np.ndarray:
        '''
        Evaluates all criteria related to a single ply.
        
        Parameters
        ------------------
        ply: PlyProperty
        
            Ply material properties.
            
        plyStresses: ndarray [nSCply]
        
            Ply stresses in local reference frame
        
        oldPlyIndexes: ndarray [NFI]
        
            History failure indexes of this ply

        limitFIDen: bool
        
            Flag for denominator limitation for damage propagation
        
        Returns
        ------------------
        plyIndexes: ndarray [NFI]
        
            Failure indexes
        
        '''
        stressVector = np.zeros(6)
        plyIndexes   = np.zeros_like(oldPlyIndexes)

        #* Assemble stress vector and plane angles
        if self.isElement3D:
            
            stressVector = plyStresses.copy()
            
            matrixAngles = np.array([ 0, 5, 10, 15, 30, 45, ply.a0_degree, 60, 75, 90, 
                                        105, 120, 180-ply.a0_degree, 135, 150, 165, 170]) / 180.0*np.pi

            kinkAngles = np.linspace(0, np.pi, int(180/self.KINK_ANGLE_SPACING)+1, endpoint=True)
            
        else:
            
            stressVector[self.n11] = plyStresses[self.n11_2d]
            stressVector[self.n22] = plyStresses[self.n22_2d]
            stressVector[self.n12] = plyStresses[self.n12_2d]

            matrixAngles = np.array([0.0, ply.a0])
            
            kinkAngles = np.array([0.0, 0.5*np.pi])
    
        #* ------------------------------------------------------
        #* 1) Matrix cracking
        plyIndexes[0] = self.get_matrix_cracking(ply, stressVector, oldPlyIndexes[0], matrixAngles, limitFIDen)
        
        #* 3) Fibre tension
        plyIndexes[2] = self.get_fibre_tension(ply, stressVector, oldPlyIndexes[2])
            
        #* 2) Matrix splitting and 4) Fibre Kinking
        plyIndexes[1], plyIndexes[3] = self.get_matrix_splitting_n_fibre_kinking(ply, stressVector, 
                                                    oldPlyIndexes[1], oldPlyIndexes[3], kinkAngles, limitFIDen)
        
        #* 5) NCF criteria: matrix interface @ 22.11.2019
        if self.isElement3D and ply.Zt is not None:
            
            plyIndexes[4] = self.get_matrix_interface(ply, stressVector, oldPlyIndexes[4])
        
        return plyIndexes
    

if __name__ == '__main__': 

    stress_vector = [1.724570E-02, -2.794602E-01, -7.775985E-02, -8.645547E-02, -1.544686E-02, -1.859897E-02]
    true_uvarm = [1.652488E-03,  0.000000E+00,  6.736600E-06,  0.000000E+00,  0.000000E+00,  1.652488E-03,  1.001000E+03]
    old_uvarm = np.zeros(7)

    uvarm_3d = UVARM(nSCply=6, material='IM7/8551-7')

    cur_uvarm = uvarm_3d.get_uvarm(stress_vector, old_uvarm, limitFIDen=False)
    
    print(cur_uvarm)
    print(cur_uvarm-np.array(true_uvarm))
    