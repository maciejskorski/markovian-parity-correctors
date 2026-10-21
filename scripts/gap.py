"""Gap at n=3,4: Zhou-Bruck's claimed exponent vs the truth. (a) observed binary chains, (b) hidden chain (A,C,..,C,A)."""
import itertools
import numpy as np

print("(a) binary chain, next-bit bias <= beta; grid over (d1,d2) in the allowed box, both x0")
print(f"{'n':>2} {'beta':>5} {'max|bias(S_n)|':>15} {'beta^ceil(n/2)':>15} {'ZB claim b2^((n-1)/2)':>22}")
for n, beta in itertools.product((3, 4), (0.2, 0.4, 0.6, 0.8)):
    g = np.linspace((1 - beta) / 2, (1 + beta) / 2, 81)
    sup = 0
    for d1, d2 in itertools.product(g, g):
        K = np.array([[d1, 1 - d1], [1 - d2, d2]]) * [1, -1]              # K = P diag(1,-1)
        sup = max(sup, np.abs(np.linalg.matrix_power(K, n) @ [1, 1]).max())  # |(K^n 1)_{x0}|
    b2 = 2 * beta / (1 + beta**2)                                          # two-sided bias valid for ZB's hypothesis
    assert abs(sup - beta ** ((n + 1) // 2)) < 1e-9
    print(f"{n:2d} {beta:5.1f} {sup:15.6f} {beta**((n+1)//2):15.6f} {b2**((n-1)/2):22.6f}")
print("=> truth = beta^ceil(n/2) (attained at d1=d2=(1+-beta)/2): ZB's bound is valid here but loose\n")

print("(b) hidden chain: theta_1=A, theta_j=(A,C_j), theta_n=A; bits x_1=x_n=A, x_j=C_j; A,C_j fair independent")
for n in (3, 4):
    paths = [(a, *c, a) for a in (0, 1) for c in itertools.product((0, 1), repeat=n - 2)]
    th = lambda x, j: (x[0],) if j in (1, n) else (x[0], x[j - 1])      # hidden state at time j (1-based)
    e = 0
    for i, j, k in itertools.combinations(range(1, n + 1), 3):
        for key in {(th(x, i), th(x, k)) for x in paths}:
            sel = [x[j - 1] for x in paths if (th(x, i), th(x, k)) == key]
            e = max(e, 2 * abs(np.mean(sel) - .5))                       # |P[x_j=1|theta_i,theta_k]-1/2| <= e/2
    bias = np.mean([(-1) ** (x[0] ^ x[-1]) for x in paths])
    print(f"n={n}: noise hypothesis holds with e = {e}; yet bias(x_1^x_n) = {bias} (ZB lemma, q=2: <= e^(1/2) = {e**.5})")
    assert e == 0 and bias == 1
