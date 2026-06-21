#!/usr/bin/env python3
"""Diagnostic script to trace alpha calculation."""
import mpmath as mp

def main():
    mp.mp.dps = 100
    print("mpmath version:", mp.__version__)

    pi = mp.pi
    e = mp.e
    tau = mp.log(pi / e)
    gamma_val = 196 * tau / (49 + tau)
    alpha0 = (pi - e)**2 / (pi**2 * mp.sqrt(2*pi))
    eps0 = mp.mpf(1) / (28**2)

    print("pi =", pi)
    print("e =", e)
    print("tau =", tau)
    print("gamma =", gamma_val)
    print("alpha0 =", alpha0)
    print("eps0 =", eps0)

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
            if g != 0:
                val *= mp.exp(-eps0)
            T[i,j] = val
        return T

    a = mp.mpf(1)/137
    for step in range(5):  # just a few steps for diagnostics
        T = build_T(a)
        inv = (mp.eye(9)-T)**-1
        Z_M = inv[1,1]
        anew = alpha0 / (abs(Z_M)**2)
        print(f"step {step}: a = {a}, Z_M = {Z_M}, anew = {anew}")
        a = anew

    print("Final alpha:", float(a))
    print("1/alpha =", float(1/a))

if __name__ == "__main__":
    main()
