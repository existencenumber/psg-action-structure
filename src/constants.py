"""
Fundamental transcendental numbers and derived constants.
All numbers are computed at the precision set by mpmath.mp.dps.
"""
import mpmath as mp

def init_constants():
    """Initialize and return fundamental constants."""
    pi = mp.pi
    e = mp.e
    tau = mp.log(pi / e)
    gamma = 196 * tau / (49 + tau)          # Euler constant from torsion equation
    alpha0 = (pi - e)**2 / (pi**2 * mp.sqrt(2*pi))  # tree-level fine-structure
    eps0 = mp.mpf(1) / (28**2)              # regularization parameter
    return pi, e, tau, gamma, alpha0, eps0
