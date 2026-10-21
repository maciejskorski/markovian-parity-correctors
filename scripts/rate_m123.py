"""Experiment: max of the per-bit rate rho(K) over order-m binary chains with next-bit bias |b_s| <= beta (m=1,2,3),
vs the conjectured sharp rate beta^(1/(m+1)) and the nilpotent bounds r_* (root of r^L = beta(1+..+r^(L-1))), (L beta)^(1/L)."""
import itertools
import numpy as np

def kernel(b, m):                                     # K(s,s') = ((-1)^x + b_s)/2 for s' = (s_2..s_m, x)
    S = list(itertools.product((0, 1), repeat=m)); ix = {s: i for i, s in enumerate(S)}; K = np.zeros((2**m, 2**m))
    for s in S:
        for x in (0, 1): K[ix[s], ix[s[1:] + (x,)]] = ((-1)**x + b[ix[s]]) / 2
    return K

rho = lambda b, m: np.abs(np.linalg.eigvals(kernel(b, m))).max()
rng = np.random.default_rng(0)

def r_star(beta, L):
    lo, hi = 0, 1
    for _ in range(80):
        r = (lo + hi) / 2
        lo, hi = (r, hi) if r**L < beta * sum(r**k for k in range(L)) else (lo, r)
    return lo

print(f"{'m':>2} {'beta':>6} {'max rho':>9} {'beta^(1/L)':>11} {'ratio':>6} {'r_*':>7} {'(L beta)^(1/L)':>15}  argmax (b_s/beta)")
for m in (1, 2, 3):
    L, n = m + 1, 2**m
    for beta in (0.01, 0.05, 0.1, 0.3, 0.5, 0.8):
        best, arg = max(((rho(beta * np.array(v), m), beta * np.array(v)) for v in itertools.product((-1, 1), repeat=n)), key=lambda t: t[0])
        vbest = best                                   # best over vertices b_s = +-beta
        for _ in range(60):                            # random-restart hill climbing over the whole box
            b = rng.uniform(-beta, beta, n); r = rho(b, m); step = beta / 2
            for _ in range(400):
                c = np.clip(b + step * rng.standard_normal(n), -beta, beta); rc = rho(c, m)
                if rc > r: b, r = c, rc
                else: step *= .985
            if r > best: best, arg = r, b
        tag = "vertex" if best <= vbest + 1e-12 else "INTERIOR"
        print(f"{m:2d} {beta:6.2f} {best:9.5f} {beta**(1/L):11.5f} {best/beta**(1/L):6.3f} {r_star(beta, L):7.4f} {min(1, (L*beta)**(1/L)):15.4f}  {tag} {np.round(arg/beta, 2)}")

# threshold beta_m where the conjecture rho <= beta^(1/L) starts to fail (vertex maximizers, bisection)
ratio = lambda beta, m: max(rho(beta * np.array(v), m) for v in itertools.product((-1, 1), repeat=2**m)) / beta**(1 / (m + 1))
print()
for m in (1, 2, 3):
    if ratio(0.99, m) <= 1 + 1e-9: print(f"m={m}: ratio = 1 for all tested beta (conjecture holds up to beta=0.99)"); continue
    lo, hi = 1e-4, 0.99
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if ratio(mid, m) <= 1 + 1e-9 else (lo, mid)
    print(f"m={m}: conjecture rho <= beta^(1/(m+1)) holds for beta <= {lo:.4f} and fails right after (ratio {ratio(hi, m):.4f})")
