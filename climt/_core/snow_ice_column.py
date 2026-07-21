"""Reusable 1-D implicit snow/ice/soil heat-diffusion column solver.

Index convention: node 0 is the BOTTOM (ice base / soil interface), node
n-1 is the TOP (atmosphere-facing surface). ``rho``/``c``/``kappa`` are
length n-1 arrays defined on the layers between nodes.
"""
from dataclasses import dataclass
import numpy as np

from .backend import jit_compile
from .tridiagonal import solve_tridiagonal, tridiagonal_matvec


@dataclass
class Dirichlet:
    value: float


@dataclass
class Flux:
    value: float  # W/m^2, positive downward INTO the column at that boundary


# Integer boundary-condition tags for the jit kernel (dataclasses / isinstance
# cannot cross the njit boundary).
_BC_DIRICHLET = 0
_BC_FLUX = 1


@jit_compile
def _solve_column_kernel(rho, c, kappa, temperature, dt, dz,
                         top_flag, top_val, bot_flag, bot_val):
    """Jitable core of :func:`solve_column`.

    Assembles the Crank-Nicolson tridiagonal system with the two boundary
    rows folded into the diagonals, then solves it with the shared Thomas
    solver. Mathematically identical to the previous ``scipy.sparse``/
    ``spsolve`` implementation.
    """
    num_layers = temperature.shape[0]
    K_mid = kappa
    K_interface = 0.5 * (kappa[:-1] + kappa[1:])
    heat_capacity = rho * c
    heat_capacity_int = 0.5 * (heat_capacity[:-1] + heat_capacity[1:])
    mu_inv_int = dt / (heat_capacity_int * 2.0 * dz * dz)

    r = np.zeros(num_layers)
    a_sub = np.zeros(num_layers)
    a_sup = np.zeros(num_layers)
    r[1:-1] = K_interface * mu_inv_int
    dp = 1.0 + 2.0 * r
    dm = 1.0 - 2.0 * r
    a_sub[:-2] = -mu_inv_int * K_mid[:-1]
    a_sup[2:] = -mu_inv_int * K_mid[1:]

    # RHS = mat_rhs @ temperature, where mat_rhs has diagonal dm and
    # off-diagonals -a_sub / -a_sup. Matching scipy.spdiags' offset
    # convention (verified: sub-diagonal = a_sub[:-1], super-diagonal =
    # a_sup[1:]), the banded layout used by solve_tridiagonal takes the
    # sub-diagonal from -a_sub[:-1] and the super-diagonal from -a_sup[1:].
    rhs = tridiagonal_matvec(-a_sub[:-1], dm, -a_sup[1:], temperature)

    dp_lhs = dp.copy()
    a_sub_lhs = a_sub.copy()
    a_sup_lhs = a_sup.copy()

    # Top boundary (node n-1)
    if top_flag == _BC_DIRICHLET:
        dp_lhs[-1] = 1.0
        a_sub_lhs[-2] = 0.0
        rhs[-1] = top_val
    else:  # Flux (Neumann): K dT/dz = -flux at the top face
        dp_lhs[-1] = -1.0
        a_sub_lhs[-2] = 1.0
        rhs[-1] = -top_val * dz / K_mid[-1]

    # Bottom boundary (node 0)
    if bot_flag == _BC_DIRICHLET:
        dp_lhs[0] = 1.0
        a_sup_lhs[1] = 0.0
        rhs[0] = bot_val
    else:  # Flux into the base
        dp_lhs[0] = -1.0
        a_sup_lhs[1] = 1.0
        rhs[0] = -bot_val * dz / K_mid[0]

    return solve_tridiagonal(a_sub_lhs[:-1], dp_lhs, a_sup_lhs[1:], rhs)


def solve_column(rho, c, kappa, temperature, dt, dz, top_bc, bottom_bc):
    """Solve one implicit snow/ice/soil heat-diffusion column.

    ``top_bc``/``bottom_bc`` are :class:`Dirichlet` or :class:`Flux`
    instances; they are translated to integer tags + scalar values before
    calling the jitable :func:`_solve_column_kernel`.
    """
    top_flag = _BC_DIRICHLET if isinstance(top_bc, Dirichlet) else _BC_FLUX
    bot_flag = _BC_DIRICHLET if isinstance(bottom_bc, Dirichlet) else _BC_FLUX
    return _solve_column_kernel(
        np.asarray(rho, dtype=float),
        np.asarray(c, dtype=float),
        np.asarray(kappa, dtype=float),
        np.asarray(temperature, dtype=float),
        float(dt),
        float(dz),
        top_flag,
        float(top_bc.value),
        bot_flag,
        float(bottom_bc.value),
    )
