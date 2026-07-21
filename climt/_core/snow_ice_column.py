"""Reusable 1-D implicit snow/ice/soil heat-diffusion column solver.

Index convention: node 0 is the BOTTOM (ice base / soil interface), node
n-1 is the TOP (atmosphere-facing surface). ``rho``/``c``/``kappa`` are
length n-1 arrays defined on the layers between nodes.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


@dataclass
class Dirichlet:
    value: float


@dataclass
class Flux:
    value: float  # W/m^2, positive downward INTO the column at that boundary


def solve_column(rho, c, kappa, temperature, dt, dz, top_bc, bottom_bc):
    num_layers = temperature.shape[0]
    K_mid = kappa
    K_interface = 0.5 * (kappa[:-1] + kappa[1:])
    heat_capacity = rho * c
    heat_capacity_int = 0.5 * (heat_capacity[:-1] + heat_capacity[1:])
    mu_inv_int = dt / (heat_capacity_int * 2 * dz * dz)

    r = np.zeros(num_layers)
    a_sub = np.zeros(num_layers)
    a_sup = np.zeros(num_layers)
    r[1:-1] = K_interface * mu_inv_int
    dp = 1 + 2 * r
    dm = 1 - 2 * r
    a_sub[:-2] = -mu_inv_int * K_mid[:-1]
    a_sup[2:] = -mu_inv_int * K_mid[1:]

    mat_lhs = sparse.spdiags([a_sub, dp, a_sup], [-1, 0, 1],
                             num_layers, num_layers, format="csc")
    mat_rhs = sparse.spdiags([-a_sub, dm, -a_sup], [-1, 0, 1],
                             num_layers, num_layers, format="csc")
    rhs = mat_rhs * temperature

    # Top boundary (node n-1)
    if isinstance(top_bc, Dirichlet):
        mat_lhs[-1, -1] = 1; mat_lhs[-1, -2] = 0
        rhs[-1] = top_bc.value
    else:  # Flux (Neumann): K dT/dz = -flux at the top face
        mat_lhs[-1, -1] = -1; mat_lhs[-1, -2] = 1
        rhs[-1] = -top_bc.value * dz / K_mid[-1]

    # Bottom boundary (node 0)
    if isinstance(bottom_bc, Dirichlet):
        mat_lhs[0, 0] = 1; mat_lhs[0, 1] = 0
        rhs[0] = bottom_bc.value
    else:  # Flux into the base
        mat_lhs[0, 0] = -1; mat_lhs[0, 1] = 1
        rhs[0] = -bottom_bc.value * dz / K_mid[0]

    return spsolve(mat_lhs, rhs)
