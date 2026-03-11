# RVE modeling of open hole plate

This document summarizes the RVE modeling of open hole plate with different hole sizes.

Boundary conditions are applied via constraint equations to homogenize the
three-dimensional representative volume elements (RVE) to an orthotropic material.

The effective 6x6 stiffness matrix $C$ (Voigt notation) is as follows,
where $1, 2, 3$ correspond to $x, y, z$ directions, respectively.

- $\boldsymbol{\sigma} = C\boldsymbol{\epsilon}$
- $\boldsymbol{\sigma} = [\sigma_{11}, \sigma_{22}, \sigma_{33}, \sigma_{23}, \sigma_{13}, \sigma_{12}]^T$
- $\boldsymbol{\epsilon} = [\epsilon_{11}, \epsilon_{22}, \epsilon_{33}, \gamma_{23}, \gamma_{13}, \gamma_{12}]^T$
- $S = C^{-1}$ (compliance matrix)
- $\gamma_{ij} = 2 \cdot \epsilon_{ij}$ (engineering shear strain, $i \neq j$)

## 1. Boundary conditions

Denote the dimensions of the RVE as $l_{x}$, $l_{y}$, $l_{z}$, and the hole radius as $r_\text{hole}$.

Denote the representative points as $\text{RP}$ to implement the boundary conditions.

Denote the displacement as $u_1, u_2, u_3$ for $x, y, z$ directions, respectively.

Use representative points to introduce the strain components, i.e.,
$$
\begin{align}
\epsilon_{11} &= u_1(\text{RP}_{11}), \\
\epsilon_{22} &= u_1(\text{RP}_{22}), \\
\epsilon_{33} &= u_1(\text{RP}_{33}), \\
\gamma_{23} &= u_1(\text{RP}_{23}), \\
\gamma_{13} &= u_1(\text{RP}_{13}), \\
\gamma_{12} &= u_1(\text{RP}_{12}).
\end{align}
$$

### 1.1 Periodic boundary condition (PBC)

The periodic boundary conditions are implemented to the pairs of nodes in the opposite faces of the RVE,
which are called the master face and the slave face.

Denote the master face in $x, y, z$ directions as $\text{Mx}, \text{My}, \text{Mz}$,
and the slave face as $\text{Sx}, \text{Sy}, \text{Sz}$.

Assuming small deformation, the periodic boundary conditions are:

$$
\begin{align}
& u_1(\text{Mx}) - u_1(\text{Sx}) =& l_x \cdot \epsilon_{11} \\
& u_2(\text{Mx}) - u_2(\text{Sx}) =& 0.5 \cdot l_x \cdot \gamma_{12} \\
& u_3(\text{Mx}) - u_3(\text{Sx}) =& 0.5 \cdot l_x \cdot \gamma_{13}
\end{align}
$$

$$
\begin{align}
& u_1(\text{My}) - u_1(\text{Sy}) =& 0.5 \cdot l_y \cdot \gamma_{12} \\
& u_2(\text{My}) - u_2(\text{Sy}) =& l_y \cdot \epsilon_{22} \\
& u_3(\text{My}) - u_3(\text{Sy}) =& 0.5 \cdot l_y \cdot \gamma_{23}
\end{align}
$$

$$
\begin{align}
& u_1(\text{Mz}) - u_1(\text{Sz}) =& 0.5 \cdot l_z \cdot \gamma_{13} \\
& u_2(\text{Mz}) - u_2(\text{Sz}) =& 0.5 \cdot l_z \cdot \gamma_{23} \\
& u_3(\text{Mz}) - u_3(\text{Sz}) =& l_z \cdot \epsilon_{33}
\end{align}
$$

### 1.2 Linear boundary condition (LBC)

The linear boundary conditions force the displacement of nodes on RVE faces to be linear distributions.

Assuming small deformation, one corner of the RVE is at the origin ($x=0, y=0, z=0$),
and the faces goes through the origin are called the slave faces,
the boundary conditions are:

$$
\begin{align}
& u_1(\text{Mx}) =& \epsilon_{11} \cdot l_x + 0.5 \cdot (\gamma_{12} \cdot y(\text{Mx}) + \gamma_{13} \cdot z(\text{Mx})) \\
& u_2(\text{Mx}) =& \epsilon_{22} \cdot y(\text{Mx}) + 0.5 \cdot (\gamma_{12} \cdot l_x + \gamma_{23} \cdot z(\text{Mx})) \\
& u_3(\text{Mx}) =& \epsilon_{33} \cdot z(\text{Mx}) + 0.5 \cdot (\gamma_{13} \cdot l_x + \gamma_{23} \cdot y(\text{Mx}))
\end{align}
$$

$$
\begin{align}
& u_1(\text{My}) =& \epsilon_{11} \cdot x(\text{My}) + 0.5 \cdot (\gamma_{12} \cdot l_y + \gamma_{13} \cdot z(\text{My})) \\
& u_2(\text{My}) =& \epsilon_{22} \cdot l_y + 0.5 \cdot (\gamma_{12} \cdot x(\text{My}) + \gamma_{23} \cdot z(\text{My})) \\
& u_3(\text{My}) =& \epsilon_{33} \cdot z(\text{My}) + 0.5 \cdot (\gamma_{23} \cdot l_y + \gamma_{13} \cdot x(\text{My}))
\end{align}
$$

$$
\begin{align}
& u_1(\text{Mz}) =& \epsilon_{11} \cdot x(\text{Mz}) + 0.5 \cdot (\gamma_{13} \cdot l_z + \gamma_{12} \cdot y(\text{Mz})) \\
& u_2(\text{Mz}) =& \epsilon_{22} \cdot y(\text{Mz}) + 0.5 \cdot (\gamma_{23} \cdot l_z + \gamma_{12} \cdot x(\text{Mz})) \\
& u_3(\text{Mz}) =& \epsilon_{33} \cdot l_z + 0.5 \cdot (\gamma_{13} \cdot x(\text{Mz}) + \gamma_{23} \cdot y(\text{Mz}))
\end{align}
$$

$$
\begin{align}
& u_1(\text{Sx}) =& 0.5 \cdot (\gamma_{12} \cdot y(\text{Sx}) + \gamma_{13} \cdot z(\text{Sx})) \\
& u_2(\text{Sx}) =& \epsilon_{22} \cdot y(\text{Sx}) + 0.5 \cdot (\gamma_{23} \cdot z(\text{Sx})) \\
& u_3(\text{Sx}) =& \epsilon_{33} \cdot z(\text{Sx}) + 0.5 \cdot (\gamma_{23} \cdot y(\text{Sx}))
\end{align}
$$

$$
\begin{align}
& u_1(\text{Sy}) =& \epsilon_{11} \cdot x(\text{Sy}) + 0.5 \cdot (\gamma_{13} \cdot z(\text{Sy})) \\
& u_2(\text{Sy}) =& 0.5 \cdot (\gamma_{12} \cdot x(\text{Sy}) + \gamma_{23} \cdot z(\text{Sy})) \\
& u_3(\text{Sy}) =& \epsilon_{33} \cdot z(\text{Sy}) + 0.5 \cdot (\gamma_{13} \cdot x(\text{Sy}))
\end{align}
$$

$$
\begin{align}
& u_1(\text{Sz}) =& \epsilon_{11} \cdot x(\text{Sz}) + 0.5 \cdot (\gamma_{12} \cdot y(\text{Sz})) \\
& u_2(\text{Sz}) =& \epsilon_{22} \cdot y(\text{Sz}) + 0.5 \cdot (\gamma_{12} \cdot x(\text{Sz})) \\
& u_3(\text{Sz}) =& 0.5 \cdot (\gamma_{13} \cdot x(\text{Sz}) + \gamma_{23} \cdot y(\text{Sz}))
\end{align}
$$

### 1.3 Linear boundary condition with Z-direction periodicity (LBC-2)

The LBC-2 force the displacement of nodes on $x, y$ direction faces to be linear distributions,
and the displacement of nodes on $z$ direction faces is periodic.

## 2. Homogenized material properties

### 2.1 Raw material properties

| Material   | $E_{11}$ | $E_{22}$ | $E_{33}$ | $G_{23}$ | $G_{13}$ | $G_{12}$ | $\nu_{23}$ | $\nu_{13}$ | $\nu_{12}$ |
| ---------- | -------- | -------- | -------- | -------- | -------- | -------- | ---------- | ---------- | ---------- |
| steel      | 210000   | 210000   | 210000   | 80770    | 80770    | 80770    | 0.3        | 0.3        | 0.3        |
| IM7/8551-7 | 165000   | 8400     | 8400     | 2800     | 5600     | 5600     | 0.5        | 0.34       | 0.34       |

### 2.2 Big hole

The geometric parameters are:

- $l_{x}=l_{y}=30$ mm
- $l_z=10$ mm
- $r_\text{hole}=10$ mm

| Case  | $E_{11}$ | $E_{22}$ | $E_{33}$ | $G_{23}$ | $G_{13}$ | $G_{12}$ | $\nu_{23}$ | $\nu_{13}$ | $\nu_{12}$ |
| ----- | -------- | -------- | -------- | -------- | -------- | -------- | ---------- | ---------- | ---------- |
| steel | 210000   | 210000   | 210000   | 80770    | 80770    | 80770    | 0.3        | 0.3        | 0.3        |
| PBC   | 94675    | 94675    | 136732   | 38886    | 38886    | 20078    | 0.209      | 0.209      | 0.189      |
| LBC-2 | 96790    | 96790    | 136738   | 39977    | 39977    | 37573    | 0.213      | 0.213      | 0.200      |
| LBC   | 110613   | 110613   | 136796   | 48619    | 48619    | 43599    | 0.242      | 0.242      | 0.246      |

| Case       | $E_{11}$ | $E_{22}$ | $E_{33}$ | $G_{23}$ | $G_{13}$ | $G_{12}$ | $\nu_{23}$ | $\nu_{13}$ | $\nu_{12}$ |
| ---------- | -------- | -------- | -------- | -------- | -------- | -------- | ---------- | ---------- | ---------- |
| IM7/8551-7 | 165000   | 8400     | 8400     | 2800     | 5600     | 5600     | 0.5        | 0.34       | 0.34       |
| PBC        | 96853    | 6385     | 8407     | 2158     | 3948     | 2118     | 0.379      | 0.199      | 0.220      |
| LBC-2      | 97389    | 6433     | 8407     | 2284     | 3987     | 4354     | 0.382      | 0.200      | 0.243      |

### 2.3 Middle-sized hole

The geometric parameters are:

- $l_{x}=l_{y}=100$ mm
- $l_z=10$ mm
- $r_\text{hole}=10$ mm

| Case  | $E_{11}$ | $E_{22}$ | $E_{33}$ | $G_{23}$ | $G_{13}$ | $G_{12}$ | $\nu_{23}$ | $\nu_{13}$ | $\nu_{12}$ |
| ----- | -------- | -------- | -------- | -------- | -------- | -------- | ---------- | ---------- | ---------- |
| steel | 210000   | 210000   | 210000   | 80770    | 80770    | 80770    | 0.3        | 0.3        | 0.3        |
| PBC   | 189618   | 189618   | 203215   | 75856    | 75856    | 73737    | 0.282      | 0.282      | 0.299      |
| LBC-2 | 189947   | 189947   | 203222   | 75870    | 75870    | 74332    | 0.282      | 0.282      | 0.299      |
| LBC   | 200680   | 200680   | 203410   | 77873    | 77873    | 77337    | 0.295      | 0.295      | 0.297      |

| Case       | $E_{11}$ | $E_{22}$ | $E_{33}$ | $G_{23}$ | $G_{13}$ | $G_{12}$ | $\nu_{23}$ | $\nu_{13}$ | $\nu_{12}$ |
| ---------- | -------- | -------- | -------- | -------- | -------- | -------- | ---------- | ---------- | ---------- |
| IM7/8551-7 | 165000   | 8400     | 8400     | 2800     | 5600     | 5600     | 0.5        | 0.34       | 0.34       |
| PBC        | 145897   | 8119     | 8400     | 2738     | 5367     | 5348     | 0.483      | 0.300      | 0.354      |
| LBC-2      | 146005   | 8120     | 8400     | 2740     | 5367     | 5419     | 0.483      | 0.300      | 0.354      |

### 2.4 Small hole

The geometric parameters are:

- $l_{x}=l_{y}=400$ mm
- $l_z=10$ mm
- $r_\text{hole}=10$ mm

| Case  | $E_{11}$ | $E_{22}$ | $E_{33}$ | $G_{23}$ | $G_{13}$ | $G_{12}$ | $\nu_{23}$ | $\nu_{13}$ | $\nu_{12}$ |
| ----- | -------- | -------- | -------- | -------- | -------- | -------- | ---------- | ---------- | ---------- |
| steel | 210000   | 210000   | 210000   | 80770    | 80770    | 80770    | 0.3        | 0.3        | 0.3        |
| PBC   | 208498   | 208498   | 209568   | 80453    | 80453    | 80326    | 0.298      | 0.298      | 0.300      |
| LBC-2 | 208501   | 208501   | 209568   | 80453    | 80453    | 80328    | 0.298      | 0.298      | 0.300      |
| LBC   | 209417   | 209417   | 209588   | 80588    | 80588    | 80554    | 0.299      | 0.299      | 0.299      |

| Case       | $E_{11}$ | $E_{22}$ | $E_{33}$ | $G_{23}$ | $G_{13}$ | $G_{12}$ | $\nu_{23}$ | $\nu_{13}$ | $\nu_{12}$ |
| ---------- | -------- | -------- | -------- | -------- | -------- | -------- | ---------- | ---------- | ---------- |
| IM7/8551-7 | 165000   | 8400     | 8400     | 2800     | 5600     | 5600     | 0.5        | 0.34       | 0.34       |
| PBC        | 163158   | 8380     | 8400     | 2796     | 5584     | 5586     | 0.498      | 0.336      | 0.342      |
| LBC-2      | 163159   | 8380     | 8400     | 2796     | 5584     | 5586     | 0.498      | 0.336      | 0.342      |

## 3. Summary

The following conclusions are drawn:

