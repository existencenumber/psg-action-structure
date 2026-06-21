"""
Compute all physical constants from the self-consistent alpha.
Reference: Main text, Table I; Supplemental Material, Sec. IV.
"""
import mpmath as mp
from .constants import init_constants
from .transfer_matrix import build_T

def compute_all(alpha=None, verbose=True):
    """
    Compute all core physical constants.
    
    Parameters:
        alpha: if None, solve self-consistently first
        verbose: print results
    
    Returns:
        dict of constants
    """
    pi, e, tau, gamma, alpha0, eps0 = init_constants()
    
    if alpha is None:
        from .self_consistent import solve_alpha
        alpha = solve_alpha(tol=mp.mpf('1e-30'), verbose=verbose)
    
    T = build_T(alpha, eps0)
    inv = (mp.eye(9) - T) ** (-1)
    
    # Closed-path correction factors
    W_A = inv[0, 0] - 1   # Additive
    W_M = inv[1, 1] - 1   # Multiplicative
    W_B = inv[6, 6] - 1   # Braiding
    
    # Fine-structure constant (self-consistency check)
    alpha_check = alpha0 / (abs(1 + W_M) ** 2)
    
    # Weinberg angle
    C_ratio = mp.mpf(1) / (alpha ** (-mp.mpf('0.5')))  # C_A / C_M = sqrt(alpha)
    tan2_W = (abs(1 + W_M) ** 2) / (abs(1 + W_A) ** 2) * (C_ratio ** 2) / 3
    sin2_W = tan2_W / (1 + tan2_W)
    
    # Strong coupling
    pi_minus_e = pi - e
    alpha_s0 = alpha * pi**2 * mp.exp(pi_minus_e) * (1 + pi_minus_e / (pi + e))
    alpha_s = alpha_s0 / (abs(1 + W_B) ** 2)
    
    # Muon-to-electron mass ratio
    phi_r = 3 * pi * alpha / 2
    phi_ang = mp.mpf(1) / 2
    phi_tot = mp.sqrt(phi_r**2 + phi_ang**2)
    m_ratio = (1 - mp.cos(phi_tot)) / (1 - mp.cos(phi_r))
    
    # Hierarchy (Higgs VEV / Planck mass)
    hierarchy = mp.exp(-gamma / (2 * alpha))
    
    # Dark energy density
    M_Pl_eV = mp.mpf('1.22e28')  # Planck mass in eV
    N_eff = 14 * (1 + gamma / (4 * pi))
    rho_L = M_Pl_eV**4 * (gamma / (2 * pi))**2 * (e / pi) ** (N_eff / alpha)
    
    # Inflation
    n_s = 1 - mp.mpf(3) / 55          # N_* = 55
    r = mp.mpf(2) / (55 ** mp.mpf('1.5'))
    
    results = {
        'alpha': alpha,
        'alpha_inv': 1 / alpha,
        'alpha_s': alpha_s,
        'sin2_W': sin2_W,
        'm_ratio': m_ratio,
        'hierarchy': hierarchy,
        'rho_L': rho_L,
        'n_s': n_s,
        'r': r,
    }
    
    if verbose:
        print("=" * 60)
        print("Process Space Geometry - Physical Constants")
        print("=" * 60)
        print(f"alpha^(-1)        = {float(1/alpha):.12f}")
        print(f"alpha_s(MZ)       = {float(alpha_s):.6f}")
        print(f"sin^2(theta_W)   = {float(sin2_W):.6f}")
        print(f"m_mu / m_e       = {float(m_ratio):.6f}")
        print(f"v / M_Pl         = {float(hierarchy):.3e}")
        print(f"rho_Lambda (eV^4) = {float(rho_L):.3e}")
        print(f"n_s              = {float(n_s):.4f}")
        print(f"r                = {float(r):.5f}")
        print("=" * 60)
    
    return results
