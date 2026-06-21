"""
Quantum capacities of the nine domains in the broken state.
Reference: Supplemental Material, Eq. (S15).
"""
import mpmath as mp

def get_capacities(alpha):
    """
    Return list of 9 quantum capacities.
    Order: A, M, I, D, S, P, B, H, C
    """
    C = [mp.mpf(0)] * 9
    C[0] = mp.mpf(1)                       # Additive
    C[1] = alpha ** (-mp.mpf('0.5'))       # Multiplicative
    C[2] = alpha ** (-mp.mpf('0.25'))      # Integral
    C[3] = alpha ** (-mp.mpf('0.25'))      # Differential
    C[4] = alpha ** (-mp.mpf('0.5'))       # Spectral (= Multiplicative)
    C[5] = alpha ** (-mp.mpf(1) / 3)       # Functional Integral
    C[6] = mp.sqrt(2)                      # Braiding (topological constant)
    C[7] = mp.mpf(1)                       # Homotopy
    C[8] = mp.mpf(1)                       # Category
    return C
