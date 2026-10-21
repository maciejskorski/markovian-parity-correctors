"""Exact proof that max |bias(S_n)| = beta^2 for n=3,4 over binary chains with next-bit bias <= beta.
Param.: d1=(1+b*u)/2, d2=(1+b*v)/2, b=beta in [0,1], u,v in [-1,1]; bias(S_n)=(K^n 1)_{x0} is divisible by b^2.
f = bias/b^2 is a polynomial; 1-f >= 0 and 1+f >= 0 on the cube is certified by nonnegative Bernstein
coefficients (exact rationals, adaptive subdivision); f = +-1 identically in beta along a vertex edge (d1=d2=(1+-beta)/2) shows it is attained."""
import itertools, sys
from fractions import Fraction as Fr
from math import comb
import sympy as sp

S = sp.symbols("s t w")                                   # unit-cube coords: b=s, u=2t-1, v=2w-1
s_, t_, w_ = S
b, u, v = s_, 2 * t_ - 1, 2 * w_ - 1
K = sp.Matrix([[(1 + b*u)/2, (1 + b*u)/2 - 1], [1 - (1 + b*v)/2, -(1 + b*v)/2]])    # K = P diag(1,-1)

def bernstein(P):
    N = [P.degree(x) for x in S]; a = {m: Fr(int(c.p), int(c.q)) for m, c in P.terms()}
    return {(i, j, k): sum(a.get((l, m, n), 0) * Fr(comb(i, l) * comb(j, m) * comb(k, n), comb(N[0], l) * comb(N[1], m) * comb(N[2], n))
                           for l in range(i + 1) for m in range(j + 1) for n in range(k + 1))
            for i in range(N[0] + 1) for j in range(N[1] + 1) for k in range(N[2] + 1)}, N

def nonneg(P, depth=10):
    c, N = bernstein(P)
    if min(c.values()) >= 0: return True
    if any(c[tuple(N[d] * e[d] for d in range(3))] < 0 for e in itertools.product((0, 1), repeat=3)): return False   # corner value < 0
    if depth == 0: raise RuntimeError("depth exhausted")
    for h in itertools.product((0, 1), repeat=3):          # split cube into 8 sub-cubes
        sub = {x: (x + hh) / 2 for x, hh in zip(S, h)}
        if not nonneg(sp.Poly(P.as_expr().subs(sub, simultaneous=True), *S), depth - 1): return False
    return True

for n in map(int, sys.argv[1:] or ["3", "4"]):
    for x0 in (0, 1):
        f = sp.cancel(sp.expand((K**n * sp.ones(2, 1))[x0]) / b**2)
        assert f.is_polynomial(*S)
        ok = all(nonneg(sp.Poly(sp.expand(1 + sg * f), *S)) for sg in (1, -1))
        top = [(tt, ww) for tt, ww in itertools.product((0, 1), repeat=2) if sp.expand(f.subs({t_: tt, w_: ww})) in (1, -1)]
        print(f"n={n} x0={x0}: |bias(S_n)| <= beta^2 certified: {ok}; equality for all beta at (t,w) = {top[:2]}")
        assert ok and top
print("=> max |bias(S_n)| = beta^2 for n in", sys.argv[1:] or [3, 4])
