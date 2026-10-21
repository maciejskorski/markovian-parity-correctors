"""XOR bias trajectories n -> bias(S_n): twisted-kernel formula vs Monte Carlo (and the Cayley-Hamilton recursion)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N, nmax, x0, rng = 10**6, 14, 0, np.random.default_rng(1)
configs = {"monotone $(\\delta_1,\\delta_2)=(0.8,0.7)$": (0.8, 0.7),
           "oscillating $(\\delta_1,\\delta_2)=(0.2,0.3)$": (0.2, 0.3)}
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
for ax, (name, (d1, d2)) in zip(axes, configs.items()):
    K = np.array([[d1, d1 - 1], [1 - d2, -d2]])               # K = P Z
    tau, lam = d1 - d2, d1 + d2 - 1
    assert np.allclose(K @ K, tau * K + lam * np.eye(2))      # Cayley-Hamilton
    v, exact = np.ones(2), [1.0]
    for _ in range(nmax):
        v = K @ v; exact.append(v[x0])                        # bias(S_n) = (K^n 1)_{x0}
    rec = [1.0, exact[1]]
    for n in range(2, nmax + 1): rec.append(tau * rec[-1] + lam * rec[-2])
    assert np.allclose(rec, exact)
    x, par, mc = np.full(N, x0), np.zeros(N, int), [1.0]      # simulate N chains, running XOR
    for n in range(nmax):
        x = np.where(rng.random(N) < np.where(x == 0, d1, d2), x, 1 - x)
        par ^= x; mc.append(1 - 2 * par.mean())
    ax.plot(exact, "-", c="#b2182b", label="kernel $(K^n\mathbf{1})_{x_0}$")
    ax.plot(mc, "o", c="#2166ac", ms=4, label=f"Monte Carlo ($N=10^6$)")
    ax.axhline(0, c="gray", lw=.5); ax.set_title(name, fontsize=10); ax.set_xlabel("$n$")
    print(f"{name}: max |kernel - MC| = {np.max(np.abs(np.array(mc) - exact)):.1e}  (MC noise ~ {1/np.sqrt(N):.1e})")
axes[0].set_ylabel("bias$(S_n)$"); axes[0].legend(frameon=False)
fig.tight_layout(); fig.savefig("figures/xor_trajectories.png", dpi=150)
