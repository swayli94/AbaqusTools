'''
Implicit modelling of open hole.
'''
import os
import numpy as np
from typing import List, Tuple, Dict


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
    
    if dx0 <= 0 or dx1 <= 0 or dy0 <= 0 or dy1 <= 0:
        print('dx0 = %.2f, dx1 = %.2f, dy0 = %.2f, dy1 = %.2f' % (dx0, dx1, dy0, dy1))
        raise ValueError('Hole is out of bounds of the plate. Please check the geometry parameters.')
    
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


def update_parameters(parameters: dict, target_volume_fraction: float=0.4,
                    use_implicit_modelling: bool=True) -> dict:
    '''
    Update the parameters for implicit modelling and mesh-independent fastener.
    
    Note: only one fastener in the specimen.
    '''
    len_x = parameters['pGeo']['len_x_plate']
    len_y = parameters['pGeo']['len_y_plate']
    
    pFastener = parameters['pGeo']['fasteners'][0]
    r_hole = pFastener['r_hole']
    xc_hole = pFastener['x_center']
    yc_hole = pFastener['y_center']

    width_partition, vf_hole = calculate_partition_dimensions(
        len_x, len_y, r_hole, xc_hole, yc_hole, 
        target_volume_fraction=target_volume_fraction
    )
    
    if use_implicit_modelling:
    
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
        
    else:
        
        material = {
            'E11': parameters['pMesh']['E11'],
            'E22': parameters['pMesh']['E22'], 
            'E33': parameters['pMesh']['E33'],
            'G12': parameters['pMesh']['G12'],
            'G13': parameters['pMesh']['G13'],
            'G23': parameters['pMesh']['G23'],
            'nu12': parameters['pMesh']['nu12'],
            'nu13': parameters['pMesh']['nu13'],
            'nu23': parameters['pMesh']['nu23']
        }
    
    # Update the parameters for implicit modelling
    parameters['pMesh']['ImplicitModelling'] = {}
    parameters['pMesh']['ImplicitModelling']['width_partition'] = width_partition
    parameters['pMesh']['ImplicitModelling']['vf_hole'] = vf_hole
    parameters['pMesh']['ImplicitModelling'].update(material)
    
    return parameters


'''
Analytical stress field (Lekhnitskii solution) for laminate plate with a circular hole.
'''
from typing import List, Tuple
from lamkit.analysis.material import Ply, Material
from lamkit.analysis.laminate import Laminate
from lamkit.lekhnitskii.utils import generate_meshgrid
from lamkit.utils import evaluate_combined_load_plate, NUMERIC_KEYS


def read_ctf_from_rf_dat(fname: str, fastener_index: int = 0) -> Tuple[float, float, float]:
    '''
    Read connector transfer forces CTF1, CTF2, and CTF3 from a RF dat file.

    Parameters
    ----------
    fname : str
        Path to the RF dat file (e.g., 'Job_MIF_0_0-RF.dat').
    fastener_index : int
        Fastener index (default 0, matching FASTENER0_CTF1/CTF2/CTF3).

    Returns
    -------
    CTF1, CTF2, CTF3 : float
        In-plane connector transfer forces in x, y, and z [N].
    '''
    key1 = 'FASTENER%d_CTF1' % fastener_index
    key2 = 'FASTENER%d_CTF2' % fastener_index
    key3 = 'FASTENER%d_CTF3' % fastener_index
    CTF1 = CTF2 = CTF3 = None
    with open(fname, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                k, v = parts[0], float(parts[1])
                if k == key1:
                    CTF1 = v
                elif k == key2:
                    CTF2 = v
                elif k == key3:
                    CTF3 = v
    if CTF1 is None or CTF2 is None or CTF3 is None:
        raise ValueError('CTF not found in %s for fastener %d' % (fname, fastener_index))
    return CTF1, CTF2, CTF3


def calculate_bypass_load(
        X: np.ndarray, Y: np.ndarray,
        N11: np.ndarray, N22: np.ndarray, N12: np.ndarray,
        bearing_force: float, angle_degree: float,
        ) -> dict:
    '''
    Bypass load calculation for a mesh-independent fastener (HyperSizer method).

    Parameters
    ----------
    X, Y : 1-D array_like
        Element centroid coordinates [mm].  Duplicate rows are removed
        automatically (extraction artefact from extract_mid_plane_strain).
    N11, N22, N12 : 1-D array_like
        Mid-plane section forces [N/mm] at element centroids.
    bearing_force : float
        Bearing force magnitude P [N].
    angle_degree : float
        Bearing force angle in degrees, 0 means along 1 (x) direction,
        positive counter-clockwise.

    Returns
    -------
    results : dict
        N_x_bypass  : float [N/mm]  bypass section force in x
        N_y_bypass  : float [N/mm]  bypass section force in y
        N_xy_bypass : float [N/mm]  bypass section shear force
        F_x_bypass  : float [N]     bypass force resultant in x  (= N_x_bypass  * W_Yspan)
        F_y_bypass  : float [N]     bypass force resultant in y  (= N_y_bypass  * W_Xspan)
        N_x_avg     : float [N/mm]  area-averaged N11 over partition square
        N_y_avg     : float [N/mm]  area-averaged N22 over partition square
        N_xy_avg    : float [N/mm]  area-averaged N12 over partition square
        W_Xspan     : float [mm]    effective partition width  in x
        W_Yspan     : float [mm]    effective partition height in y
    '''
    X   = np.asarray(X,   dtype=float)
    Y   = np.asarray(Y,   dtype=float)
    N11 = np.asarray(N11, dtype=float)
    N22 = np.asarray(N22, dtype=float)
    N12 = np.asarray(N12, dtype=float)

    # deduplicate on (X, Y)
    ux = np.unique(X)
    uy = np.unique(Y)
    dx = float(np.mean(np.diff(ux)))
    dy = float(np.mean(np.diff(uy)))

    # effective span = centroid range + one element size
    W_Xspan = float(ux[-1] - ux[0]) + dx
    W_Yspan = float(uy[-1] - uy[0]) + dy

    # area-weighted average — uniform mesh, so simple mean
    N_x_avg  = float(np.mean(N11))
    N_y_avg  = float(np.mean(N22))
    N_xy_avg = float(np.mean(N12))

    # bearing force
    p     = bearing_force
    alpha = np.radians(angle_degree)

    # ------------------------------------------------------------------
    # HyperSizer area-average bypass formula
    #
    # The element-average section force over a partition window spanning the
    # fastener zone equals the mean of the upstream (incoming) and downstream
    # (bypass) section forces:
    #     N_avg = (N_incoming + N_bypass) / 2
    #
    # Combining this with the in-plane force balance
    #     N_incoming * w = N_bypass * w + F_bearing
    # gives:
    #     N_x_bypass  = sign(N_x_avg) * ( |N_x_avg| - P*cos(alpha) / (2*W_Yspan) )
    #     N_y_bypass  = sign(N_y_avg) * ( |N_y_avg| - P*sin(alpha) / (2*W_Xspan) )
    #     N_xy_bypass = N_xy_avg                    # shear has no bearing term
    #
    #   P, alpha   bearing force magnitude and load angle (from connector CTF)
    #   W_Yspan    partition height, used with N_x  [mm]
    #   W_Xspan    partition width,  used with N_y  [mm]
    #   sign(.)    preserves direction for compression cases (N_avg < 0)
    #
    # Self-check: by construction, 2 * N_x_avg * W_Yspan - 2 * N_x_bypass * W_Yspan = CTF1
    # ------------------------------------------------------------------

    N_x_bypass_hs  = float(np.sign(N_x_avg) * (abs(N_x_avg) - p * np.cos(alpha) / (2.0 * W_Yspan)))
    N_y_bypass_hs  = float(np.sign(N_y_avg) * (abs(N_y_avg) - p * np.sin(alpha) / (2.0 * W_Xspan)))
    N_xy_bypass_hs = N_xy_avg

    return {
        'N_x_bypass'    : N_x_bypass_hs,
        'N_y_bypass'    : N_y_bypass_hs,
        'N_xy_bypass'   : N_xy_bypass_hs,
        'F_x_bypass'    : N_x_bypass_hs  * W_Yspan,
        'F_y_bypass'    : N_y_bypass_hs  * W_Xspan,
        'N_x_avg'       : N_x_avg,
        'N_y_avg'       : N_y_avg,
        'N_xy_avg'      : N_xy_avg,

        # --- geometry ---
        'W_Xspan' : W_Xspan,
        'W_Yspan' : W_Yspan,
    }


def calculate_ply_level_stress_field(pMesh: dict,
                    r_hole: float, total_thickness: float,
                    N11: float, N22: float, N12: float,
                    bearing_load: float = 0.0, angle_load_degree: float = 0.0,
                    characteristic_distance: float = 1.0,
                    n_points_radial: int = 5, n_points_angular: int = 32,
                    ) -> Tuple[dict, List[dict]]:
    '''
    Calculate the ply-level stress field around the open hole using the Lekhnitskii solution.
    
    Parameters
    ----------
    pMesh : dict
        The mesh parameters for the problem.
    r_hole : float
        The radius of the hole.
    total_thickness : float
        The total thickness of the laminate.
    N11, N22, N12 : float
        The applied in-plane forces per unit length at infinity, in the 11, 22 and 12 directions, respectively.
    bearing_load : float, default 0.0
        The bearing load applied to the hole.
    angle_load_degree : float, default 0.0
        The angle of the applied load with respect to the 1-2 coordinate system (degrees).
    characteristic_distance : float, default 1.0
        The characteristic distance for failure evaluation, used for mesh generation.
    n_points_radial : int, default 5
        The number of points in the radial direction for the mesh grid.
    n_points_angular : int, default 32
        The number of points in the angular direction for the mesh grid.
    calculate_failure : bool, default False
        Whether to calculate the LaRC05 failure indices and modes for each ply.
        
    Returns
    -------
    field : dict
        The stress and strain fields at the mid-plane of the laminate,
        with keys 'X', 'Y', 'sigma_xx', 'sigma_yy', 'tau_xy', 'epsilon_x', 'epsilon_y', 'gamma_xy'.
    results_by_plies : List[dict]
        The stress and strain fields for each ply,
        with keys 'sigma_1', 'sigma_2', 'tau_12'.
    '''
    pIM = pMesh['ImplicitModelling']
    properties = {
        'E11': pIM['E11'], 'E22': pIM['E22'],
        'nu12': pIM['nu12'], 'G12': pIM['G12']
    }
    layup = pMesh['plate_CompositePly_orientationValue']
    if pMesh['plate_CompositeLayup_symmetric']:
        layup = layup + layup[::-1]
    ply_thickness = pMesh['composite_ply_thickness']
    if ply_thickness * len(layup) != total_thickness:
        raise ValueError(f'Total laminate thickness {ply_thickness * len(layup)} exceeds plate thickness {total_thickness}.')
    
    material = Material(name='implicit_modelling_material',
                    properties=properties, check_larc05=False)
    ply = Ply(material=material, thickness=ply_thickness)
    laminate = Laminate(stacking=layup, plies=ply)

    mesh = generate_meshgrid(
        hole_radius=r_hole,
        plate_radius=r_hole+characteristic_distance,
        n_points_radial=n_points_radial,
        n_points_angular=n_points_angular,
        radial_cluster_power=2.0,
        )

    results_by_plies_2, mid_plane_field = evaluate_combined_load_plate(
        laminate=laminate,
        sigma_xx_inf=N11/total_thickness,
        sigma_yy_inf=N22/total_thickness,
        tau_xy_inf=N12/total_thickness,
        load=bearing_load,
        angle_load_degree=angle_load_degree,
        hole_radius=r_hole,
        thickness=total_thickness,
        x=mesh["X"], y=mesh["Y"]
        )
    
    # Average the top and bottom surface fields of each ply
    results_by_plies = [{} for _ in range(laminate.n_ply)]
    for index_ply in range(laminate.n_ply):
        result_b = results_by_plies_2[2*index_ply]
        result_t = results_by_plies_2[2*index_ply + 1]
        for key in NUMERIC_KEYS:
            results_by_plies[index_ply][key] = 0.5*(result_b[key] + result_t[key])
    
    field = {
        "X": mesh["X"],
        "Y": mesh["Y"],
        "z_edges": np.array(laminate.z_position, dtype=float) + 0.5*total_thickness,
        "sigma_xx": mid_plane_field["sigma_x"],
        "sigma_yy": mid_plane_field["sigma_y"],
        "tau_xy": mid_plane_field["tau_xy"],
        "epsilon_x": mid_plane_field["epsilon_x"],
        "epsilon_y": mid_plane_field["epsilon_y"],
        "gamma_xy": mid_plane_field["gamma_xy"],
    }
    return field, results_by_plies


def read_tecplot(filename: str) -> List[Dict[str, np.ndarray]]:
    '''
    Read a point-data Tecplot ASCII file and return one dict per zone.

    Supported header formats
    ------------------------
    Variables= V1 V2 V3 ...
    zone T="name"              # name may contain spaces
    zone T="name" I=n          # optional I= (number of rows) and other params

    Each returned dict contains:
        'name' : str            zone title (T="..." stripped)
        'V1'   : np.ndarray     one key per variable, shape (n_rows,)
        'V2'   : ...

    Rows that cannot be converted to floats (e.g., extra zone headers inside
    a zone) are silently skipped.

    Parameters
    ----------
    filename : str
        Path to the Tecplot ASCII file.

    Returns
    -------
    List[Dict[str, np.ndarray]]
        One dict per zone in order of appearance.
    '''
    if not os.path.isfile(filename):
        raise FileNotFoundError('Tecplot file not found: %s' % filename)

    with open(filename, 'r') as f:
        lines = f.readlines()

    # --- parse variable names from the first non-blank line ---
    var_names: List[str] = []
    start_idx = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            continue
        low = stripped.lower()
        if low.startswith('variables'):
            # "Variables= V1 V2 ..." or "Variables = V1 V2 ..."
            after_eq = stripped.split('=', 1)[1]
            var_names = after_eq.split()
            start_idx = i + 1
            break

    if not var_names:
        raise ValueError('No "Variables=" header found in %s' % filename)

    # --- parse zones and data rows ---
    zones: List[Dict[str, np.ndarray]] = []
    current_name: str = ''
    current_rows: List[List[float]] = []

    def _flush_zone():
        if current_name and current_rows:
            arr = np.array(current_rows, dtype=float)   # (n_rows, n_vars)
            entry: Dict[str, np.ndarray] = {'name': current_name}
            for j, v in enumerate(var_names):
                entry[v] = arr[:, j]
            zones.append(entry)

    for line in lines[start_idx:]:
        stripped = line.strip()
        if not stripped:
            continue
        low = stripped.lower()

        if low.startswith('zone'):
            _flush_zone()
            current_rows = []
            # extract zone title from T="..." or T=' ...'
            import re
            m = re.search(r'[Tt]\s*=\s*["\']([^"\']*)["\']', stripped)
            current_name = m.group(1).strip() if m else stripped
        else:
            try:
                current_rows.append([float(v) for v in stripped.split()])
            except ValueError:
                pass   # skip unparseable lines

    _flush_zone()   # flush the last zone

    return zones


def save_tecplot(filename: str, field: dict, results_by_plies: List[dict]):
    '''
    Save the field and ply-level results to a Tecplot file for visualization.
    '''
    ni = field['X'].shape[0]
    nj = field['X'].shape[1]
    nk = len(results_by_plies)
    z_edges = field["z_edges"]

    with open(filename, 'w') as f:
        f.write('Variables= X Y Z S11 S22 S12\n')
        f.write('zone T=" PLATE PARTITION_CIRCLE " I=%d J=%d K=%d\n'%(ni, nj, nk))
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    f.write(' %14.6E'%(field['X'][i, j]))
                    f.write(' %14.6E'%(field['Y'][i, j]))
                    f.write(' %14.6E'%(0.5*(z_edges[k] + z_edges[k+1])))
                    f.write(' %14.6E'%(results_by_plies[k]['sigma_1'][i, j]))
                    f.write(' %14.6E'%(results_by_plies[k]['sigma_2'][i, j]))
                    f.write(' %14.6E'%(results_by_plies[k]['tau_12'][i, j]))
                    f.write('\n')

