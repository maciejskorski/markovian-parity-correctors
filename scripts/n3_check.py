"""n=3 check: bias(X1^X2^X3) from the twisted kernel vs brute-force enumeration of all 8 paths."""
import itertools
import numpy as np

rng = np.random.default_rng(0)
print(f"{'d1':>6} {'d2':>6} {'x0':>3} {'enumeration':>12} {'pi K^3 1':>12}")
for _ in range(5):
    d1, d2 = rng.random(2); x0 = int(rng.integers(2))
    P = np.array([[d1, 1 - d1], [1 - d2, d2]])        # holding probabilities
    K = P @ np.diag([1, -1])                           # twisted kernel K = P Z, Z = diag((-1)^x)
    kernel = (np.linalg.matrix_power(K, 3) @ np.ones(2))[x0]
    brute = sum(np.prod([P[a, b] for a, b in zip((x0,) + x, x)]) * (-1) ** sum(x)
                for x in itertools.product((0, 1), repeat=3))
    assert abs(kernel - brute) < 1e-12
    print(f"{d1:6.3f} {d2:6.3f} {x0:3d} {brute:12.8f} {kernel:12.8f}")
print("OK: enumeration == kernel formula")
