# climt/_components/cork/sw/kernels.py
#
# Meador & Weaver (1980) two-stream SW solver following the RTE/RRTMGP
# implementation (rte/kernels/mo_rte_solver_kernels.F90).
#
# Key references:
#   - Meador & Weaver 1980, doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
#   - Zdunkowski et al. 1980, Contributions to Atmospheric Physics 53, 147-66 (PIFM)
#   - Shonk & Hogan 2008, doi:10.1175/2007JCLI1940.1 (adding method)
import numpy as np

from ..common import njit, prange

_MIN_K = 1.0e-12
_MIN_MU0 = 1.0e-8


@njit
def _delta_scale(tau, ssa, g):
    """Apply delta-Eddington scaling to optical properties.

    Removes the forward scattering peak for better two-stream accuracy.

    Args:
        tau, ssa, g: scalars — optical depth, single scattering albedo, asymmetry

    Returns:
        tau_s, ssa_s, g_s: delta-scaled values
    """
    f = g * g
    tau_s = tau * (1.0 - ssa * f)
    ssa_s = ssa * (1.0 - f) / (1.0 - ssa * f) if (1.0 - ssa * f) > 1e-30 else 0.0
    g_s = (g - f) / (1.0 - f) if (1.0 - f) > 1e-30 else 0.0
    return tau_s, ssa_s, g_s


@njit
def _sw_dif_and_source(tau, w0, g, mu0):
    """Compute per-layer diffuse R/T and direct-beam source terms.

    Following RRTMGP's sw_dif_and_source (mo_rte_solver_kernels.F90).
    Uses Zdunkowski PIFM gamma coefficients.

    Args:
        tau: layer optical depth (scalar, delta-scaled)
        w0: single scattering albedo (scalar, delta-scaled)
        g: asymmetry parameter (scalar, delta-scaled)
        mu0: cosine of solar zenith angle (scalar)

    Returns:
        Rdif: diffuse reflectance
        Tdif: diffuse transmittance
        Rdir: fraction of direct beam scattered upward (source_up per unit flux)
        Tdir: fraction of direct beam scattered downward (source_dn per unit flux)
        Tnoscat: direct beam transmittance (Beer's law)
    """
    # PIFM gamma coefficients (Zdunkowski et al. 1980)
    gamma1 = (8.0 - w0 * (5.0 + 3.0 * g)) * 0.25
    gamma2 = 3.0 * (w0 * (1.0 - g)) * 0.25

    # k = sqrt(gamma1^2 - gamma2^2), floored to avoid div by 0 at conservative limit
    k = np.sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), _MIN_K))
    exp_minusktau = np.exp(-tau * k)
    exp_minus2ktau = exp_minusktau * exp_minusktau

    # RT_term: common factor (refactored for numerical stability)
    RT_term = 1.0 / (k * (1.0 + exp_minus2ktau) + gamma1 * (1.0 - exp_minus2ktau))

    # Diffuse reflectance and transmittance (M&W Eqs. 25, 26)
    Rdif = RT_term * gamma2 * (1.0 - exp_minus2ktau)
    Tdif = RT_term * 2.0 * k * exp_minusktau

    # Direct beam: unscattered transmittance
    mu0_s = max(mu0, _MIN_MU0)
    Tnoscat = np.exp(-tau / mu0_s)

    # Direct beam source terms (M&W Eqs. 14, 15)
    k_mu = k * mu0_s
    denom_dir = 1.0 - k_mu * k_mu
    if abs(denom_dir) < 1e-30:
        denom_dir = 1e-30

    RT_term_dir = w0 * RT_term / denom_dir

    gamma3 = (2.0 - 3.0 * mu0_s * g) * 0.25
    gamma4 = 1.0 - gamma3
    alpha1 = gamma1 * gamma4 + gamma2 * gamma3  # Eq. 16
    alpha2 = gamma1 * gamma3 + gamma2 * gamma4  # Eq. 17

    k_gamma3 = k * gamma3
    k_gamma4 = k * gamma4

    Rdir = RT_term_dir * (
        (1.0 - k_mu) * (alpha2 + k_gamma3)
        - (1.0 + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau
        - 2.0 * (k_gamma3 - alpha2 * k_mu) * exp_minusktau * Tnoscat
    )

    Tdir = -RT_term_dir * (
        (1.0 + k_mu) * (alpha1 + k_gamma4) * Tnoscat
        - (1.0 - k_mu) * (alpha1 - k_gamma4) * exp_minus2ktau * Tnoscat
        - 2.0 * (k_gamma4 + alpha1 * k_mu) * exp_minusktau
    )

    # Energy clamping (Robin Hogan / Peter Ukkonen)
    Rdir = max(0.0, min(Rdir, 1.0 - Tnoscat))
    Tdir = max(0.0, min(Tdir, 1.0 - Tnoscat - Rdir))

    return Rdif, Tdif, Rdir, Tdir, Tnoscat


@njit
def _adding(nlev, sfc_albedo, Rdif, Tdif, src_up, src_dn, src_sfc,
            flux_dn_dir):
    """Adding method for diffuse fluxes (Shonk & Hogan 2008).

    Convention: index 0 = surface, index nlev = TOA.
    Layer k sits between interface k (bottom) and k+1 (top).

    Args:
        nlev: number of layers
        sfc_albedo: surface albedo (scalar)
        Rdif: (nlev,) diffuse reflectance per layer
        Tdif: (nlev,) diffuse transmittance per layer
        src_up: (nlev,) upward source from direct beam scattering per layer
        src_dn: (nlev,) downward source from direct beam scattering per layer
        src_sfc: surface source (direct beam reflected by surface)
        flux_dn_dir: (nlev+1,) direct beam flux at interfaces

    Returns:
        flux_up: (nlev+1,) upward diffuse flux at interfaces
        flux_dn: (nlev+1,) downward diffuse flux at interfaces
    """
    # Arrays for combined albedo and source (SH08 Eqs. 9-11)
    albedo = np.zeros(nlev + 1)  # effective diffuse albedo below each interface
    src = np.zeros(nlev + 1)     # accumulated upward source below each interface
    denom = np.zeros(nlev)

    # Surface boundary
    albedo[0] = sfc_albedo
    src[0] = src_sfc

    # Upward sweep: build combined albedo and source from surface to TOA
    for k in range(nlev):
        denom[k] = 1.0 / (1.0 - Rdif[k] * albedo[k])  # Eq. 10
        # Eq. 9: combined albedo above layer k
        albedo[k + 1] = Rdif[k] + Tdif[k] * Tdif[k] * albedo[k] * denom[k]
        # Eq. 11: combined source above layer k
        src[k + 1] = src_up[k] + Tdif[k] * denom[k] * (
            src[k] + albedo[k] * src_dn[k]
        )

    # Fluxes
    flux_up = np.zeros(nlev + 1)
    flux_dn = np.zeros(nlev + 1)

    # TOA: no incoming diffuse
    flux_dn[nlev] = 0.0
    # Eq. 12 at TOA
    flux_up[nlev] = flux_dn[nlev] * albedo[nlev] + src[nlev]

    # Downward sweep: propagate diffuse fluxes from TOA to surface
    for k in range(nlev - 1, -1, -1):
        # Eq. 13
        flux_dn[k] = (
            Tdif[k] * flux_dn[k + 1] + Rdif[k] * src[k] + src_dn[k]
        ) * denom[k]
        # Eq. 12
        flux_up[k] = flux_dn[k] * albedo[k] + src[k]

    return flux_up, flux_dn


@njit
def _sw_two_stream_core(tau, ssa, asymmetry, zenith, albedo, solar_flux, weights):
    """Multi-band, multi-g-point SW radiative transfer.

    Meador & Weaver (1980) two-stream with delta-Eddington scaling and
    adding method for inter-layer coupling. Follows RTE/RRTMGP.

    Args:
        tau: (nband, ngpt, nlev, ncol) extinction optical depth
        ssa: (nband, ngpt, nlev, ncol) single scattering albedo
        asymmetry: (nband, ngpt, nlev, ncol) asymmetry parameter
        zenith: (ncol,) solar zenith angle, radians
        albedo: (ncol,) surface albedo
        solar_flux: (nband, ngpt) TOA solar flux per g-point, W/m^2
        weights: (nband, ngpt) g-point weights

    Returns:
        up_band: (nband, nlev+1, ncol) per-band upward flux
        down_band: (nband, nlev+1, ncol) per-band downward flux
        up_broad: (nlev+1, ncol) broadband upward flux
        down_broad: (nlev+1, ncol) broadband downward flux
    """
    nband, ngpt, nlev, ncol = tau.shape

    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))

    for b in range(nband):
        for g in range(ngpt):
            w = weights[b, g]
            for i in prange(ncol):
                mu0 = np.cos(zenith[i])
                if mu0 <= 1e-4:
                    continue  # nighttime

                # Per-layer properties
                Rdif = np.zeros(nlev)
                Tdif = np.zeros(nlev)
                Rdir = np.zeros(nlev)
                Tdir = np.zeros(nlev)
                Tnoscat = np.zeros(nlev)

                # Direct beam flux at interfaces (top-down)
                flux_dn_dir = np.zeros(nlev + 1)
                flux_dn_dir[nlev] = 1.0  # unit flux at TOA

                # Source terms from direct beam scattering
                src_up = np.zeros(nlev)
                src_dn = np.zeros(nlev)

                # Process layers from TOA downward
                for k in range(nlev - 1, -1, -1):
                    tau_s, ssa_s, g_s = _delta_scale(
                        tau[b, g, k, i], ssa[b, g, k, i], asymmetry[b, g, k, i]
                    )
                    Rdif[k], Tdif[k], Rdir[k], Tdir[k], Tnoscat[k] = (
                        _sw_dif_and_source(tau_s, ssa_s, g_s, mu0)
                    )

                    # Direct beam at bottom of layer k
                    flux_dn_dir[k] = Tnoscat[k] * flux_dn_dir[k + 1]

                    # Source terms: direct beam incident on this layer
                    src_up[k] = Rdir[k] * flux_dn_dir[k + 1]
                    src_dn[k] = Tdir[k] * flux_dn_dir[k + 1]

                # Surface source: direct beam reflected by surface
                src_sfc = flux_dn_dir[0] * albedo[i]

                # Adding method for diffuse fluxes
                flux_up_dif, flux_dn_dif = _adding(
                    nlev, albedo[i], Rdif, Tdif, src_up, src_dn, src_sfc,
                    flux_dn_dir
                )

                # Accumulate: scale by solar flux * mu0 * weight
                scale = solar_flux[b, g] * mu0 * w
                for k in range(nlev + 1):
                    up_band[b, k, i] += flux_up_dif[k] * scale
                    down_band[b, k, i] += (flux_dn_dir[k] + flux_dn_dif[k]) * scale

    up_broad = np.zeros((nlev + 1, ncol))
    down_broad = np.zeros((nlev + 1, ncol))
    for b in range(nband):
        for k in range(nlev + 1):
            for i in range(ncol):
                up_broad[k, i] += up_band[b, k, i]
                down_broad[k, i] += down_band[b, k, i]

    return up_band, down_band, up_broad, down_broad


def sw_two_stream(tau, ssa, asymmetry, zenith, albedo, solar_flux, weights,
                  diagnostics_level=0):
    """Multi-band, multi-g-point SW radiative transfer.

    Meador & Weaver (1980) two-stream with delta-Eddington scaling and
    adding method for inter-layer coupling.

    Args:
        tau: (nband, ngpt, nlev, ncol) extinction optical depth
        ssa: (nband, ngpt, nlev, ncol) single scattering albedo
        asymmetry: (nband, ngpt, nlev, ncol) asymmetry parameter
        zenith: (ncol,) solar zenith angle, radians
        albedo: (ncol,) surface albedo
        solar_flux: (nband, ngpt) TOA solar flux per g-point, W/m^2
        weights: (nband, ngpt) g-point weights
        diagnostics_level: 0 (fluxes only), 1 (per-layer R/T + direct beam),
            2 (additionally delta-scaled properties + combined albedo)

    Returns:
        If diagnostics_level == 0:
            (up_band, down_band, up_broad, down_broad)
        If diagnostics_level > 0:
            (up_band, down_band, up_broad, down_broad, diagnostics_dict)
    """
    if diagnostics_level == 0:
        return _sw_two_stream_core(tau, ssa, asymmetry, zenith, albedo,
                                   solar_flux, weights)

    nband, ngpt, nlev, ncol = tau.shape

    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))

    diag_Rdif = np.zeros((nband, ngpt, nlev, ncol))
    diag_Tdif = np.zeros((nband, ngpt, nlev, ncol))
    diag_Tnoscat = np.zeros((nband, ngpt, nlev, ncol))
    diag_direct = np.zeros((nband, ngpt, nlev + 1, ncol))

    if diagnostics_level >= 2:
        diag_Rdir = np.zeros((nband, ngpt, nlev, ncol))
        diag_Tdir = np.zeros((nband, ngpt, nlev, ncol))
        diag_combined_albedo = np.zeros((nband, ngpt, nlev + 1, ncol))
        diag_tau_d = np.zeros((nband, ngpt, nlev, ncol))
        diag_ssa_d = np.zeros((nband, ngpt, nlev, ncol))
        diag_g_d = np.zeros((nband, ngpt, nlev, ncol))

    for b in range(nband):
        for g in range(ngpt):
            w = weights[b, g]
            for i in range(ncol):
                mu0 = np.cos(zenith[i])
                if mu0 <= 1e-4:
                    continue

                Rdif = np.zeros(nlev)
                Tdif = np.zeros(nlev)
                Rdir = np.zeros(nlev)
                Tdir = np.zeros(nlev)
                Tnoscat = np.zeros(nlev)
                flux_dn_dir = np.zeros(nlev + 1)
                flux_dn_dir[nlev] = 1.0
                src_up = np.zeros(nlev)
                src_dn = np.zeros(nlev)

                for k in range(nlev - 1, -1, -1):
                    tau_s, ssa_s, g_s = _delta_scale(
                        tau[b, g, k, i], ssa[b, g, k, i], asymmetry[b, g, k, i]
                    )
                    Rdif[k], Tdif[k], Rdir[k], Tdir[k], Tnoscat[k] = (
                        _sw_dif_and_source(tau_s, ssa_s, g_s, mu0)
                    )
                    flux_dn_dir[k] = Tnoscat[k] * flux_dn_dir[k + 1]
                    src_up[k] = Rdir[k] * flux_dn_dir[k + 1]
                    src_dn[k] = Tdir[k] * flux_dn_dir[k + 1]

                    if diagnostics_level >= 2:
                        diag_tau_d[b, g, k, i] = tau_s
                        diag_ssa_d[b, g, k, i] = ssa_s
                        diag_g_d[b, g, k, i] = g_s

                src_sfc = flux_dn_dir[0] * albedo[i]
                flux_up_dif, flux_dn_dif = _adding(
                    nlev, albedo[i], Rdif, Tdif, src_up, src_dn, src_sfc,
                    flux_dn_dir
                )

                diag_Rdif[b, g, :, i] = Rdif
                diag_Tdif[b, g, :, i] = Tdif
                diag_Tnoscat[b, g, :, i] = Tnoscat
                diag_direct[b, g, :, i] = flux_dn_dir * solar_flux[b, g] * mu0

                if diagnostics_level >= 2:
                    diag_Rdir[b, g, :, i] = Rdir
                    diag_Tdir[b, g, :, i] = Tdir
                    comb_alb = np.zeros(nlev + 1)
                    comb_alb[0] = albedo[i]
                    for k in range(nlev):
                        d = 1.0 - Rdif[k] * comb_alb[k]
                        if abs(d) < 1e-30:
                            d = 1e-30
                        comb_alb[k + 1] = Rdif[k] + Tdif[k] * Tdif[k] * comb_alb[k] / d
                    diag_combined_albedo[b, g, :, i] = comb_alb

                scale = solar_flux[b, g] * mu0 * w
                for k in range(nlev + 1):
                    up_band[b, k, i] += flux_up_dif[k] * scale
                    down_band[b, k, i] += (flux_dn_dir[k] + flux_dn_dif[k]) * scale

    up_broad = np.zeros((nlev + 1, ncol))
    down_broad = np.zeros((nlev + 1, ncol))
    for b in range(nband):
        for k in range(nlev + 1):
            for i in range(ncol):
                up_broad[k, i] += up_band[b, k, i]
                down_broad[k, i] += down_band[b, k, i]

    diag = {
        "Rdif": diag_Rdif,
        "Tdif": diag_Tdif,
        "Tnoscat": diag_Tnoscat,
        "direct_beam_flux": diag_direct,
    }
    if diagnostics_level >= 2:
        diag["Rdir"] = diag_Rdir
        diag["Tdir"] = diag_Tdir
        diag["combined_albedo"] = diag_combined_albedo
        diag["tau_delta"] = diag_tau_d
        diag["ssa_delta"] = diag_ssa_d
        diag["g_delta"] = diag_g_d

    return up_band, down_band, up_broad, down_broad, diag
