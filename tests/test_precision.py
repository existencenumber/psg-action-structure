"""
Verify numerical stability at different precisions.

Usage:
    python tests/test_precision.py
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import mpmath as mp

if __name__ == '__main__':
    for dps in [50, 100, 200]:
        mp.mp.dps = dps
        from src.physics_constants import compute_all
        results = compute_all(verbose=False)
        print(f"dps={dps}: alpha^(-1) = {float(1/results['alpha']):.12f}")
