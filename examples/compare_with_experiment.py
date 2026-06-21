"""
Compare theoretical predictions with experimental values (Table I data).

Usage:
    python examples/compare_with_experiment.py
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import mpmath as mp
mp.dps = 100

from src.physics_constants import compute_all

# Experimental values from PDG 2022, CODATA 2021, Planck 2018
experiment = {
    'alpha_inv': 137.035999084,
    'alpha_inv_err': 2.1e-8,
    'alpha_s': 0.1180,
    'alpha_s_err': 0.0009,
    'sin2_W': 0.23122,
    'sin2_W_err': 0.00003,
    'm_ratio': 206.7682830,
    'm_ratio_err': 4.6e-6,
    'n_s': 0.9649,
    'n_s_err': 0.0042,
    'r': 0.036,  # upper bound (95% CL)
}

if __name__ == '__main__':
    results = compute_all(verbose=False)
    
    print("Constant\tTheory\t\tExperiment\tDeviation")
    print("-" * 70)
    
    # alpha^{-1}
    th = float(results['alpha_inv'])
    ex = experiment['alpha_inv']
    dev = (th - ex) / ex
    print(f"alpha^-1\t{th:.12f}\t{ex:.12f}\t{dev:.2e}")
    
    # alpha_s
    th = float(results['alpha_s'])
    ex = experiment['alpha_s']
    dev = (th - ex) / ex
    print(f"alpha_s\t\t{th:.6f}\t{ex:.6f}\t{dev:.2e}")
    
    # sin^2 theta_W
    th = float(results['sin2_W'])
    ex = experiment['sin2_W']
    dev = (th - ex) / ex
    print(f"sin^2_W\t\t{th:.6f}\t{ex:.6f}\t{dev:.2e}")
    
    # m_mu / m_e
    th = float(results['m_ratio'])
    ex = experiment['m_ratio']
    dev = (th - ex) / ex
    print(f"m_mu/m_e\t{th:.6f}\t{ex:.6f}\t{dev:.2e}")
    
    # n_s
    th = float(results['n_s'])
    ex = experiment['n_s']
    dev = (th - ex) / ex
    print(f"n_s\t\t{th:.4f}\t{ex:.4f}\t{dev:.2e}")
    
    # r
    th = float(results['r'])
    print(f"r\t\t{th:.5f}\t< {experiment['r']} (95% CL)")
