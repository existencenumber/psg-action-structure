"""
Self-consistent solution for the fine-structure constant alpha.
Reference: Eq. (4) in main text.
"""
import mpmath as mp
from .transfer_matrix import build_T
from .constants import init_constants

def solve_alpha(tol=mp.mpf('1e-30'), max_iter=30, verbose=False):
    """
    Fixed-point iteration to find the self-consistent alpha.
    
    Returns:
        alpha (mp.mpf)
    """
    _, _, _, _, alpha0, eps0 = init_constants()
    alpha = mp.mpf(1) / 137  # initial guess
    
    for n in range(max_iter):
        T = build_T(alpha, eps0)
        inv = (mp.eye(9) - T) ** (-1)
        Z_M = inv[1, 1]                      # Multiplicative domain, vertex index 1
        alpha_new = alpha0 / (abs(Z_M) ** 2)
        
        if verbose:
            print(f"Iteration {n}: alpha^(-1) = {float(1/alpha_new):.12f}")
        
        if abs(alpha_new - alpha) < tol:
            if verbose:
                print(f"Converged after {n+1} iterations.")
            return alpha_new
        
        alpha = alpha_new
    
    if verbose:
        print("Warning: maximum iterations reached.")
    return alpha
