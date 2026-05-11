'''
Implicit modelling of open hole.
'''
import numpy as np
from typing import List, Tuple


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


'''
Analytical stress field (Lekhnitskii solution) for laminate plate with a circular hole.
'''
from typing import List, Dict, Any, Tuple
from lamkit.analysis.material import Ply, Material
from lamkit.analysis.laminate import Laminate
from lamkit.lekhnitskii.utils import generate_meshgrid
from lamkit.lekhnitskii.unloaded_hole import UnloadedHole


def evaluate_unloaded_hole_plate(
        laminate: Laminate, hole_radius: float,
        sigma_xx_inf: float, sigma_yy_inf: float, tau_xy_inf: float,
        x: np.ndarray, y: np.ndarray,
        ) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    '''
    Calculate the stress field around a circular hole in an infinite elastic plate
    subjected to general two-dimensional (2-D) loading.
    
    Parameters
    ----------
    laminate: Laminate
        Laminate object (units: MPa, mm)
    sigma_xx_inf : float
        applied stress in the x-direction at infinity
    sigma_yy_inf : float
        applied stress in the y-direction at infinity
    tau_xy_inf : float
        applied shear stress at infinity
    hole_radius : float
        hole radius
    x : np.ndarray
        x locations in the cartesian coordinate system
    y : np.ndarray
        y locations in the cartesian coordinate system
        
    Returns
    -------
    results_by_plies: List[Dict[str, Any]]
        List of dictionaries, each containing the results for a ply.
    mid_plane_field: Dict[str, Any]
        Dictionary containing the results for the mid-plane.
    '''

    # Meshgrid
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    out_shape = x.shape
    x_flat = np.atleast_1d(x).ravel()
    y_flat = np.atleast_1d(y).ravel()
    n_points = x_flat.shape[0]

    # Unloaded hole solution for mid-plane strains
    solution = UnloadedHole(sigma_xx_inf, sigma_yy_inf, tau_xy_inf,
                            radius=hole_radius,
                            compliance_matrix=laminate.in_plane_compliance_matrix)
    
    mid_plane_field = solution.calculate_field_results(x_flat, y_flat, out_shape)
    
    epsilon_x = mid_plane_field['epsilon_x'].ravel() # (n_points,)
    epsilon_y = mid_plane_field['epsilon_y'].ravel() # (n_points,)
    gamma_xy = mid_plane_field['gamma_xy'].ravel() # (n_points,)
    zeros_kappa = np.zeros_like(x_flat) # (n_points,)
    epsilon0 = np.stack(
        [
            epsilon_x, epsilon_y, gamma_xy,
            zeros_kappa, zeros_kappa, zeros_kappa,
        ],
        axis=1,
    )  # (n_points, 6)

    # Ply-level field dictionary
    NUMERIC_KEYS = ['sigma_1', 'sigma_2', 'tau_12']

    # Ply-level stress field
    results_by_plies = []
    for index_ply in range(laminate.n_ply):
        results_by_plies.append({key: np.zeros(n_points) for key in NUMERIC_KEYS})

    for i in range(n_points):
        results_one_point = laminate.get_ply_level_results(epsilon0[i, :]) # [2*n_ply]
        for index_ply in range(laminate.n_ply):
            result_b = results_one_point[2*index_ply]
            result_t = results_one_point[2*index_ply + 1]
            for key in NUMERIC_KEYS:
                results_by_plies[index_ply][key][i] = 0.5*(result_b[key] + result_t[key])

    # Reshape the results to the original shape
    for ply in results_by_plies:
        for key in NUMERIC_KEYS:
            ply[key] = ply[key].reshape(out_shape)
        
    return results_by_plies, mid_plane_field

def calculate_ply_level_stress_field(parameters: dict,
                    N11: float, N22: float, N12: float,
                    characteristic_distance: float = 1.0,
                    n_points_radial: int = 5, n_points_angular: int = 32,
                    ) -> Tuple[dict, List[dict]]:
    '''
    Calculate the ply-level stress field around the open hole using the Lekhnitskii solution.
    
    Parameters
    ----------
    parameters : dict
        The input parameters for the problem, including geometry, material properties and mesh settings.
    N11, N22, N12 : float
        The applied in-plane forces per unit length at infinity, in the 11, 22 and 12 directions, respectively.
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
    length_z = parameters['pGeo']['len_z_plate']
    r_hole = parameters['pGeo']['r_hole']
    pMesh = parameters['pMesh']
    pIM = pMesh['ImplicitModelling']
    properties = {
        'E11': pIM['E11'], 'E22': pIM['E22'],
        'nu12': pIM['nu12'], 'G12': pIM['G12']
    }
    layup = pMesh['plate_CompositePly_orientationValue']
    if pMesh['plate_CompositeLayup_symmetric']:
        layup = layup + layup[::-1]
    ply_thickness = pMesh['composite_ply_thickness']
    if ply_thickness * len(layup) != length_z:
        raise ValueError(f'Total laminate thickness {ply_thickness * len(layup)} exceeds plate thickness {length_z}.')
    
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

    results_by_plies, mid_plane_field = evaluate_unloaded_hole_plate(
        laminate=laminate, hole_radius=r_hole,
        sigma_xx_inf=N11/length_z,
        sigma_yy_inf=N22/length_z,
        tau_xy_inf=N12/length_z,
        x=mesh["X"], y=mesh["Y"]
        )
    
    field = {
        "X": mesh["X"],
        "Y": mesh["Y"],
        "z_edges": np.array(laminate.z_position, dtype=float) + 0.5*length_z, # z in [0, length_z]
        "sigma_xx": mid_plane_field["sigma_x"],
        "sigma_yy": mid_plane_field["sigma_y"],
        "tau_xy": mid_plane_field["tau_xy"],
        "epsilon_x": mid_plane_field["epsilon_x"],
        "epsilon_y": mid_plane_field["epsilon_y"],
        "gamma_xy": mid_plane_field["gamma_xy"],
    }
    return field, results_by_plies

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

