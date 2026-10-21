"""Exact check that max |bias(S_n)| / beta^2 = 1 (default n=3; n=4 is too slow with sympy.solve) over all binary chains with next-bit bias <= beta.
Param.: d1=(1+b*u)/2, d2=(1+b*v)/2 with b=beta in [0,1], u,v in [-1,1]. bias(S_n) = (K^n 1)_{x0} is a polynomial,
divisible by b^2; a polynomial on a box attains its extrema at stationary points of the faces/edges/corners (exact enumeration)."""
import itertools
import sympy as sp

b, u, v = X = sp.symbols("b u v")
lo_hi = {b: (0, 1), u: (-1, 1), v: (-1, 1)}
d1, d2 = (1 + b * u) / 2, (1 + b * v) / 2
K = sp.Matrix([[d1, d1 - 1], [1 - d2, -d2]])               # K = P diag(1,-1)

def extrema(f):
    vals = []
    for pat in itertools.product((0, 1, None), repeat=3):    # each variable fixed at lo / hi / free
        free = [x for x, p in zip(X, pat) if p is None]
        sub = {x: lo_hi[x][p] for x, p in zip(X, pat) if p is not None}
        g = f.subs(sub)
        sols = sp.solve([sp.diff(g, x) for x in free], free, dict=True) if free else [{}]
        for s in sols:
            rest = [x for x in free if x not in s]                    # non-isolated stationary set: f is constant on it,
            for t in itertools.product((.25, .5, .75), repeat=len(rest)):  # so sample it (its ends lie on lower-dim strata)
                pt = {**sub, **s, **{x: lo_hi[x][0] + q * (lo_hi[x][1] - lo_hi[x][0]) for x, q in zip(rest, t)}}
                pt = {x: sp.sympify(pt[x]).subs({r: pt[r] for r in rest}) for x in X}
                if all(sp.im(sp.N(pt[x])) == 0 and lo_hi[x][0] - 1e-12 <= sp.N(pt[x]) <= lo_hi[x][1] + 1e-12 for x in X):
                    vals.append(sp.N(f.subs(pt), 30))
    return min(vals), max(vals)

import sys
for n in map(int, sys.argv[1:] or ["3"]):
    for x0 in (0, 1):
        f = sp.cancel(sp.expand((K**n * sp.ones(2, 1))[x0]) / b**2)
        assert f.is_polynomial(*X)
        lo, hi = extrema(f)
        print(f"n={n} x0={x0}: bias(S_n)/beta^2 has exact range [{sp.N(lo, 12)}, {sp.N(hi, 12)}] on the box")
        assert max(abs(lo), abs(hi)) <= 1 + sp.Float("1e-25") and max(abs(lo), abs(hi)) >= 1 - sp.Float("1e-25")
print("=> max |bias(S_n)| = beta^2 attained, for n in", list(map(int, sys.argv[1:] or ["3"])))
