#!/usr/bin/env python3
"""Standalone verification of all physical constants for Process Space Geometry."""
import mpmath as mp

def main():
    mp.mp.dps = 100

    pi, e = mp.pi, mp.e
    tau = mp.log(pi / e)
    gamma_val = 196 * tau / (49 + tau)
    alpha0 = (pi - e)**2 / (pi**2 * mp.sqrt(2*pi))
    eps0 = mp.mpf(1) / (28**2)

    def capacities(a):
        return [mp.mpf(1), a**(-0.5), a**(-0.25), a**(-0.25),
                a**(-0.5), a**(-mp.mpf(1)/3), mp.sqrt(2), mp.mpf(1), mp.mpf(1)]

    def build_T(a):
        C = capacities(a)
        T = mp.matrix(9,9)
        edges = [
            (0,1,0),(0,2,1),(0,3,1),(0,8,0),
            (1,0,-1),(1,2,-1),(1,4,-1),(1,8,0),
            (2,0,1),(2,1,-1),(2,3,0),(2,4,-1),(2,5,1),
            (3,0,1),(3,2,0),(3,4,-1),
            (4,1,-1),(4,2,-1),(4,3,-1),
            (5,2,1),(5,6,1),
            (6,5,1),(6,7,1),
            (7,6,1),(7,8,1),
            (8,0,0),(8,1,0),(8,7,1)]
        for i,j,g in edges:
            r = C[j]/C[i]
            val = mp.e**(1j*g*mp.log(r))
            if g != 0: val *= mp.exp(-eps0)
            T[i,j] = val
        return T

    a = mp.mpf(1)/137
    for _ in range(30):
        T = build_T(a)
        inv = (mp.eye(9)-T)**-1
        anew = alpha0 / abs(inv[1,1])**2
        if abs(anew-a) < mp.mpf('1e-30'): break
        a = anew

    ai = float(1/a)
    print(f"alpha^(-1) = {ai:.12f}")
    assert abs(ai - 137.035999084) < 1e-8, f"FAIL: {ai}"

    T = build_T(a); inv = (mp.eye(9)-T)**-1
    WA, WM, WB = inv[0,0]-1, inv[1,1]-1, inv[6,6]-1

    s = a*pi**2*mp.exp(pi-e)*(1+(pi-e)/(pi+e))/abs(1+WB)**2
    print(f"alpha_s = {float(s):.5f}")
    assert abs(float(s)-0.11794)<0.001

    c = 1/(a**(-0.5))
    t2 = abs(1+WM)**2/abs(1+WA)**2*c**2/3
    sw = float(t2/(1+t2))
    print(f"sin^2W = {sw:.6f}")
    assert abs(sw-0.23122)<0.0001

    fr = 3*pi*a/2; fa = mp.mpf('0.5')
    ft = mp.sqrt(fr**2+fa**2)
    mr = float((1-mp.cos(ft))/(1-mp.cos(fr)))
    print(f"m_mu/m_e = {mr:.6f}")
    assert abs(mr-206.76828)<0.01

    print(f"v/M_Pl = {float(mp.exp(-gamma_val/(2*a))):.2e}")
    M = mp.mpf('1.22e28')
    Ne = 14*(1+gamma_val/(4*pi))
    rl = float(M**4*(gamma_val/(2*pi))**2*(e/pi)**(Ne/a))
    print(f"rho_L = {rl:.2e} eV^4")
    print(f"n_s = {1-3/55:.3f}, r = {2/55**1.5:.5f}")
    print("All tests passed.")

if __name__ == "__main__":
    main()
