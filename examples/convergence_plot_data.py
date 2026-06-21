"""
Generate convergence data for fixed-point iteration (Fig. 2 data).

Usage:
    python examples/convergence_plot_data.py
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import mpmath as mp
mp.mp.dps = 100

from src.constants import init_constants
from src.transfer_matrix import build_T

def convergence_data():
    _, _, _, _, alpha0, eps0 = init_constants()
    alpha = mp.mpf(1) / 137
    data = []
    for k in range(25):
        data.append(float(1 / alpha))
        T = build_T(alpha, eps0)
        inv = (mp.eye(9) - T) ** (-1)
        Z_M = inv[1, 1]
        alpha = alpha0 / (abs(Z_M) ** 2)
    return data

if __name__ == '__main__':
    data = convergence_data()
    for i, val in enumerate(data):
        print(f"{i}\t{val:.12f}")
