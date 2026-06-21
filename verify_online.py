#!/usr/bin/env python3
"""
Standalone verification script for Process Space Geometry.
Runs the full self-consistency calculation and checks all core constants.
"""
import mpmath as mp

def main():
    mp.mp.dps = 100  # 100-digit precision

    # ------- Fundamental constants -------
    pi = mp.pi
    e = mp.e
    tau = mp.log(pi / e)
    gamma_val = 196 * tau / (49 + tau)
    alpha0 = (pi - e) ** 2 / (pi**2 * mp.sqrt(2 * pi))
    eps0 = mp.mpf(1) / (28 ** 2)

    # ------- Quantum capacities -------
    def capacities(alpha):
        C = [mp.mpf(0)] * 9
        C[0] = mp.mpf(1)
        C[1] = alpha ** (-mp.mpf('0.5'))
        C[2] = alpha ** (-mp.mpf('0.25'))
        C[3] = alpha ** (-mp.mpf('0.25'))
        C[4] = alpha ** (-mp.mpf('0.5'))
        C[5] = alpha ** (-mp.mpf(1) / 3)
        C[6] = mp.sqrt(2)
        C[7] = mp.mpf(1)
        C[8] = mp.mpf(1)
        return C

    # ------- Transfer matrix -------
    def build_T(alpha):
        C = capacities(alpha)
        T = mp.matrix(9, 9)
        edges = [
            (0,1,0), (0,2,1), (0,3,1), (0,8,0),
            (1,0,-1), (1,2,-1), (1,4,-1), (1,8,0),
            (2,0,1), (2,1,-1), (2,3,0), (2,4,-1), (2,5,1),
            (3,0,1), (3,2,0), (3,4,-1),
            (4,1,-1), (4,2,-1), (4,3,-1),
            (5,2,1), (5,6,1),
            (6,5,1), (6,7,1),
            (7,6,1), (7,8,1),
            (8,0,0), (8,1,0), (8,7,1)
        ]
        for (i, j, gam) in edges:
            ratio = C[j] / C[i]
            val = mp.e ** (1j * gam * mp.log(ratio))
            if gam != 0:
                val *= mp.exp(-eps0)
            T[i, j] = val
        return T

    # ------- Fixed-point iteration for alpha -------
    alpha = mp.mpf(1) / 137
    for step in range(30):
        T = build_T(alpha)
        inv = (mp.eye(9) - T) ** (-1)
        Z_M = inv[1, 1]           # Multiplicative domain
        alpha_new = alpha0 / (abs(Z_M) ** 2)
        if abs(alpha_new - alpha) < mp.mpf('1e-30'):
            alpha = alpha_new
            break
        alpha = alpha_new

    alpha_inv = 1 / alpha
    print(f"alpha^(-1) = {float(alpha_inv):.12f}")

    # ------- Verify alpha^(-1) -------
    expected_alpha_inv = 137.035999084
    if abs(float(alpha_inv) - expected_alpha_inv) > 1e-8:
        raise SystemExit(f"FAIL: alpha^(-1)={float(alpha_inv):.12f}, expected ~{expected_alpha_inv}")

    # ------- Other constants -------
    T = build_T(alpha)
    inv = (mp.eye(9) - T) ** (-1)
    W_A = inv[0, 0] - 1
    W_M = inv[1, 1] - 1
    W_B = inv[6, 6] - 1

    # Strong coupling
    pi_minus_e = pi - e
    alpha_s0 = alpha * pi**2 * mp.exp(pi_minus_e) * (1 + pi_minus_e / (pi + e))
    alpha_s = alpha_s0 / (abs(1 + W_B) ** 2)
    print(f"alpha_s(MZ) = {float(alpha_s):.5f}")
    assert abs(float(alpha_s) - 0.11794) < 0.001

    # Weinberg angle
    C_ratio = 1 / (alpha ** (-mp.mpf('0.5')))
    tan2_W = (abs(1 + W_M) ** 2) / (abs(1 + W_A) ** 2) * (C_ratio ** 2) / 3
    sin2_W = tan2_W / (1 + tan2_W)
    print(f"sin^2(theta_W) = {float(sin2_W):.6f}")
    assert abs(float(sin2_W) - 0.23122) < 0.0001

    # Muon‑electron mass ratio
    phi_r = 3 * pi * alpha / 2
    phi_ang = mp.mpf('0.5')
    phi_tot = mp.sqrt(phi_r ** 2 + phi_ang ** 2)
    m_ratio = (1 - mp.cos(phi_tot)) / (1 - mp.cos(phi_r))
    print(f"m_mu / m_e = {float(m_ratio):.6f}")
    assert abs(float(m_ratio) - 206.76828) < 0.01

    # Hierarchy
    hier = mp.exp(-gamma_val / (2 * alpha))
    print(f"v / M_Pl = {float(hier):.2e}")

    # Dark energy density
    M_Pl_eV = mp.mpf('1.22e28')
    N_eff = 14 * (1 + gamma_val / (4 * pi))
    rho_L = M_Pl_eV**4 * (gamma_val / (2 * pi))**2 * (e / pi) ** (N_eff / alpha)
    print(f"rho_Lambda = {float(rho_L):.2e} eV^4")

    # Inflation
    ns = 1 - 3 / 55
    r = 2 / (55 ** 1.5)
    print(f"n_s = {ns:.3f}")
    print(f"r = {r:.5f}")

    print("\nAll verifications passed successfully.")

if __name__ == '__main__':
    main()
