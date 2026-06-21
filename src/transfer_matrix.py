"""
Construction of the 9x9 transfer matrix T.
Reference: Supplemental Material, Eq. (S16)-(S17).
"""
import mpmath as mp
from .capacities import get_capacities

def build_T(alpha, eps0):
    """
    Build the 9x9 transfer matrix T.
    
    Parameters:
        alpha: fine-structure constant
        eps0:  regularization parameter (1/28^2)
    
    Returns:
        mp.matrix(9,9)
    """
    C = get_capacities(alpha)
    T = mp.matrix(9, 9)
    
    # Edges: (from, to, gamma)
    # Vertex indices: 0=A,1=M,2=I,3=D,4=S,5=P,6=B,7=H,8=C
    edges = [
        # From A (index 0)
        (0, 1, 0),   # exp, free
        (0, 2, 1),   # R, limit
        (0, 3, 1),   # Q, limit
        (0, 8, 0),   # IC, free
        # From M (index 1)
        (1, 0, -1),  # log, analytic
        (1, 2, -1),  # LD, analytic
        (1, 4, -1),  # M_e, analytic
        (1, 8, 0),   # IC', free
        # From I (index 2)
        (2, 0, 1),   # R^{-1}, limit
        (2, 1, -1),  # LD^{-1}, analytic
        (2, 3, 0),   # FTC, free
        (2, 4, -1),  # L, analytic
        (2, 5, 1),   # FP, limit
        # From D (index 3)
        (3, 0, 1),   # Q^{-1}, limit
        (3, 2, 0),   # FTC, free
        (3, 4, -1),  # F, analytic
        # From S (index 4)
        (4, 1, -1),  # M_e^{-1}, analytic
        (4, 2, -1),  # L^{-1}, analytic
        (4, 3, -1),  # F^{-1}, analytic
        # From P (index 5)
        (5, 2, 1),   # FP^{-1}, limit
        (5, 6, 1),   # BT, topological
        # From B (index 6)
        (6, 5, 1),   # BT^{-1}, topological
        (6, 7, 1),   # BH, topological
        # From H (index 7)
        (7, 6, 1),   # BH^{-1}, topological
        (7, 8, 1),   # MC, topological
        # From C (index 8)
        (8, 0, 0),   # IC^{-1}, free
        (8, 1, 0),   # IC'^{-1}, free
        (8, 7, 1),   # MC^{-1}, topological
    ]
    
    for (i, j, gam) in edges:
        ratio = C[j] / C[i]
        val = mp.e ** (1j * gam * mp.log(ratio))
        if gam != 0:
            val *= mp.exp(-eps0)
        T[i, j] = val
    
    return T
