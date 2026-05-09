'''
Implicit modelling of open hole.
'''
import numpy as np
from typing import Tuple


def calculate_partition_dimensions(len_x, len_y, r_hole,
                    xc_hole, yc_hole, target_volume_fraction=0.3,
                    ratio_square=0.6):
    '''
    Calculate the width of partition square around the open hole.
    '''
    dx0 = xc_hole - r_hole
    dx1 = len_x - xc_hole - r_hole
    dy0 = yc_hole - r_hole
    dy1 = len_y - yc_hole - r_hole
    
    # volume_fraction = np.pi*r_hole**2 / width_partition**2
    width_partition = np.sqrt(np.pi*r_hole**2 / target_volume_fraction)
        
    width_partition = 2* min(
        0.5*width_partition, 
        r_hole + ratio_square*dx0, r_hole + ratio_square*dx1,
        r_hole + ratio_square*dy0, r_hole + ratio_square*dy1
        )
    
    volume_fraction = np.pi*r_hole**2 / width_partition**2
    
    return width_partition, volume_fraction

def _assemble_stiffness_compliance_matrices(
        E11: float, E22: float, E33: float,
        G12: float, G13: float, G23: float,
        nu12: float, nu13: float, nu23: float
        ) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Voigt order: (11, 22, 33, 23, 13, 12).
    Convention: nu_ij = -S_ij / S_ii.
    '''
    S = np.zeros((6, 6))
    
    S[0, 0] = 1.0 / E11
    S[1, 1] = 1.0 / E22
    S[2, 2] = 1.0 / E33
    S[3, 3] = 1.0 / G23
    S[4, 4] = 1.0 / G13
    S[5, 5] = 1.0 / G12
    
    S[0, 1] = S[1, 0] = -nu12 / E11
    S[0, 2] = S[2, 0] = -nu13 / E11
    S[1, 2] = S[2, 1] = -nu23 / E22
    
    C = np.linalg.inv(S)
    return C, S

def _voigt_reuss_averages(C: np.ndarray, S: np.ndarray, suffix='') -> dict:
    '''
    Compute Voigt/Reuss bulk and shear moduli and anisotropy indices from
    the 6x6 stiffness (C) and compliance (S) matrices in Voigt notation.
    '''
    C11, C22, C33 = C[0, 0], C[1, 1], C[2, 2]
    C12, C13, C23 = C[0, 1], C[0, 2], C[1, 2]
    C44, C55, C66 = C[3, 3], C[4, 4], C[5, 5]

    K_V = (C11 + C22 + C33 + 2*(C12 + C13 + C23)) / 9.0
    G_V = (C11 + C22 + C33 - C12 - C13 - C23 + 3*(C44 + C55 + C66)) / 15.0

    S11, S22, S33 = S[0, 0], S[1, 1], S[2, 2]
    S12, S13, S23 = S[0, 1], S[0, 2], S[1, 2]
    S44, S55, S66 = S[3, 3], S[4, 4], S[5, 5]

    K_R = 1.0 / (S11 + S22 + S33 + 2*(S12 + S13 + S23))
    G_R = 15.0 / (4*(S11 + S22 + S33) - 4*(S12 + S13 + S23) + 3*(S44 + S55 + S66))

    A_U  = 5.0 * G_V / G_R + K_V / K_R - 6.0
    A_CB = (G_V - G_R) / (G_V + G_R)
    
    result = {
        f'K_V{suffix}': K_V,
        f'G_V{suffix}': G_V,
        f'K_R{suffix}': K_R,
        f'G_R{suffix}': G_R,
        f'A_U{suffix}': A_U,
        f'A_CB{suffix}': A_CB,
    }

    return result


class OpenHoleElement(object):
    '''
    The homogenised constitutive model for elements with an open hole in the center.
    
    The plate is in the x-y plane, and the hole is along the z-axis.
    The material is orthotropic, the 1/2/3 directions are along the x/y/z axes, respectively.
    
    Parameters
    ----------
    E11, E22, E33 : float
        Young's moduli of the material.
    G12, G13, G23 : float
        Shear moduli of the material.
    nu12, nu13, nu23 : float
        Poisson's ratios of the material.
    vf_hole : float
        Volume fraction of the hole in the element.
    '''

    def __init__(self, E11: float, E22: float, E33: float,
                G12: float, G13: float, G23: float,
                nu12: float, nu13: float, nu23: float,
                vf_hole: float):

        self.material_E11 = E11
        self.material_E22 = E22
        self.material_E33 = E33
        self.material_G12 = G12
        self.material_G13 = G13
        self.material_G23 = G23
        self.material_nu12 = nu12
        self.material_nu13 = nu13
        self.material_nu23 = nu23
        self.vf_hole = vf_hole
        
        self._material_stiffness_matrix = None
        self._material_compliance_matrix = None
        self._material_A_CB = None
        self._material_ratio_E = None

        self.E11 : float = 0.0
        self.E22 : float = 0.0
        self.E33 : float = 0.0
        self.G12 : float = 0.0
        self.G13 : float = 0.0
        self.G23 : float = 0.0
        self.nu12 : float = 0.0
        self.nu13 : float = 0.0
        self.nu23 : float = 0.0
        
        self._stiffness_matrix = None
        self._compliance_matrix = None

    @property
    def material_stiffness_matrix(self) -> np.ndarray:
        '''
        Material stiffness matrix [C], Voigt notation: (11, 22, 33, 23, 13, 12).
        '''
        if self._material_stiffness_matrix is None:
            cc, ss = _assemble_stiffness_compliance_matrices(
                self.material_E11, self.material_E22, self.material_E33,
                self.material_G12, self.material_G13, self.material_G23,
                self.material_nu12, self.material_nu13, self.material_nu23
            )
            self._material_stiffness_matrix = cc
            self._material_compliance_matrix = ss
        return self._material_stiffness_matrix

    @property
    def material_compliance_matrix(self) -> np.ndarray:
        '''
        Material compliance matrix [S], Voigt notation: (11, 22, 33, 23, 13, 12).
        '''
        if self._material_compliance_matrix is None:
            cc, ss = _assemble_stiffness_compliance_matrices(
                self.material_E11, self.material_E22, self.material_E33,
                self.material_G12, self.material_G13, self.material_G23,
                self.material_nu12, self.material_nu13, self.material_nu23
            )
            self._material_stiffness_matrix = cc
            self._material_compliance_matrix = ss
        return self._material_compliance_matrix

    @property
    def material_A_CB(self) -> float:
        '''
        Material Chung-Buessem anisotropy index, ranges in [0, 1),
        where 0 is isotropic and 1 is maximally anisotropic.
        '''
        if self._material_A_CB is None:
            result = _voigt_reuss_averages(self.material_stiffness_matrix, self.material_compliance_matrix)
            self._material_A_CB = result['A_CB']
        return self._material_A_CB

    @property
    def material_ratio_E(self) -> float:
        '''
        Material ratio of elastic moduli, max(E_ii) / min(E_ii).
        '''
        if self._material_ratio_E is None:
            self._material_ratio_E = max(self.material_E11, self.material_E22, self.material_E33) / min(self.material_E11, self.material_E22, self.material_E33)
        return self._material_ratio_E

    @property
    def stiffness_matrix(self) -> np.ndarray:
        '''
        Stiffness matrix [C], Voigt notation: (11, 22, 33, 23, 13, 12).
        '''
        if self._stiffness_matrix is None:
            self.get_homogenised_properties()
            cc, ss = _assemble_stiffness_compliance_matrices(
                self.E11, self.E22, self.E33,
                self.G12, self.G13, self.G23,
                self.nu12, self.nu13, self.nu23
            )
            self._stiffness_matrix = cc
            self._compliance_matrix = ss
        return self._stiffness_matrix
    
    @property
    def compliance_matrix(self) -> np.ndarray:
        '''
        Compliance matrix [S], Voigt notation: (11, 22, 33, 23, 13, 12).
        '''
        if self._compliance_matrix is None:
            self.get_homogenised_properties()
            cc, ss = _assemble_stiffness_compliance_matrices(
                self.E11, self.E22, self.E33,
                self.G12, self.G13, self.G23,
                self.nu12, self.nu13, self.nu23
            )
            self._stiffness_matrix = cc
            self._compliance_matrix = ss
        return self._compliance_matrix

    def _get_e11(self) -> float:
        return self.material_E11 * (np.exp(-2.2*self.vf_hole) - 7*self.material_A_CB*np.clip(self.vf_hole, 0, 0.02))
    
    def _get_e22(self) -> float:
        return self.material_E22 * np.exp(-2.2*self.vf_hole)
    
    def _get_e33(self) -> float:
        '''Accurate'''
        return self.material_E33 * (1 - self.vf_hole)
    
    def _get_e11_e22_ratio(self) -> float:
        return (self.material_E11 / self.material_E22) ** 0.95
    
    def _get_e22_e11_ratio(self) -> float:
        return self.material_E22 / self.material_E11
    
    def _get_g12(self) -> float:
        # return self.material_G12 * (np.exp(-4.85*self.vf_hole**1.28) + (0.04*self.material_A_CB-0.01)*20*np.clip(self.vf_hole, 0, 0.05))
        return self.material_G12 * (1.0004423 - np.sin((1.8618281*self.material_A_CB + 1.8618281*self.material_nu12 + 0.17850223*self.material_ratio_E**0.2668298)*np.sin(self.vf_hole)/(self.material_A_CB + self.material_nu12)))

    def _get_g13(self) -> float:
        return self.material_G13 * (1.2*0.2**self.vf_hole - 0.2)
    
    def _get_g23(self) -> float:
        return self.material_G23 * (1.2*0.2**self.vf_hole - 0.2)
    
    def _get_nu12_over_e11(self) -> float:
        return self.material_nu12 / self.material_E11 * (1 + 0.5*self.vf_hole)
    
    def _get_nu13_over_e11(self) -> float:
        '''Accurate'''
        return self.material_nu13 / self.material_E11 * (6.18**self.vf_hole - self.vf_hole)
    
    def _get_nu23_over_e22(self) -> float:
        '''Accurate'''
        return self.material_nu23 / self.material_E22 / (1.0 - self.vf_hole)
    
    def _get_nu12_times_e22(self) -> float:
        return self.material_nu12 * self.material_E22 * (0.019**self.vf_hole)
        # return self.material_nu12 * self.material_E22 * (np.exp(-2.6675065*self.vf_hole*(1.0*self.nu12 + np.exp(np.sin(self.material_A_CB - 0.013364636*self.material_ratio_E)))))
    
    def _get_nu12_directly(self) -> float:
        return self.material_nu12 * 0.17296413**self.vf_hole
    
    def _get_nu12(self) -> float:
        # return self._get_nu12_over_e11() * self._get_e11() # Large error
        return self._get_nu12_times_e22() / self._get_e22()
        # return self._get_nu12_directly()
    
    def _get_nu13(self) -> float:
        return self._get_nu13_over_e11() * self._get_e11()
    
    def _get_nu23(self) -> float:
        return self._get_nu23_over_e22() * self._get_e22()

    def get_homogenised_properties(self):
        '''
        Get the homogenised material properties of the open hole element.
        The homogenised properties are calculated by a simple rule of mixtures, 
        which is not accurate but can be used as a first approximation.
        '''

        self.E11 = self._get_e11()
        self.E22 = self._get_e22()
        self.E33 = self._get_e33()
        self.G12 = self._get_g12()
        self.G13 = self._get_g13()
        self.G23 = self._get_g23()
        self.nu12 = self._get_nu12()
        self.nu13 = self._get_nu13()
        self.nu23 = self._get_nu23()

        return {
            'E11': self.E11,
            'E22': self.E22, 
            'E33': self.E33,
            'G12': self.G12,
            'G13': self.G13,
            'G23': self.G23,
            'nu12': self.nu12,
            'nu13': self.nu13,
            'nu23': self.nu23
        }


def update_parameters(parameters: dict, target_volume_fraction: float=0.4) -> dict:
    '''
    Update the parameters for implicit modelling.
    '''
    len_x = parameters['pGeo']['len_x_plate']
    len_y = parameters['pGeo']['len_y_plate']
    r_hole = parameters['pGeo']['r_hole']
    xc_hole = parameters['pGeo']['xr_hole_center'] * len_x
    yc_hole = parameters['pGeo']['yr_hole_center'] * len_y
    
    width_partition, vf_hole = calculate_partition_dimensions(
        len_x, len_y, r_hole, xc_hole, yc_hole, 
        target_volume_fraction=target_volume_fraction
    )
    
    open_hole_element = OpenHoleElement(
        E11=parameters['pMesh']['E11'],
        E22=parameters['pMesh']['E22'],
        E33=parameters['pMesh']['E33'],
        G12=parameters['pMesh']['G12'],
        G13=parameters['pMesh']['G13'],
        G23=parameters['pMesh']['G23'],
        nu12=parameters['pMesh']['nu12'],
        nu13=parameters['pMesh']['nu13'],
        nu23=parameters['pMesh']['nu23'],
        vf_hole=vf_hole
    )
    
    material = open_hole_element.get_homogenised_properties()
    
    # Update the parameters for implicit modelling
    parameters['pMesh']['ImplicitModelling'] = {}
    parameters['pMesh']['ImplicitModelling']['width_partition'] = width_partition
    parameters['pMesh']['ImplicitModelling']['vf_hole'] = vf_hole
    parameters['pMesh']['ImplicitModelling'].update(material)
    
    return parameters
