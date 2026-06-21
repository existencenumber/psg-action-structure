"""
Example: compute all physical constants and compare with experiment.

Usage:
    python examples/compute_all_constants.py
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import mpmath as mp
mp.mp.dps = 100

from src.physics_constants import compute_all

if __name__ == '__main__':
    results = compute_all(verbose=True)
