# python plot.py   -> compares godunov.csv (numerical) with the exact solution computed by your riemann_at
import numpy as np
import matplotlib.pyplot as plt

d = np.genfromtxt("godunov.csv", delimiter=",", names=True)
g = 1.4
e    = d["p"]       / ((g - 1) * d["rho"])
e_ex = d["p_exact"] / ((g - 1) * d["rho_exact"])

fig, ax = plt.subplots(2, 2, figsize=(9, 6.5), sharex=True)

panels = [("density", d["rho"], d["rho_exact"]), ("velocity", d["u"], d["u_exact"]),
          ("pressure", d["p"], d["p_exact"]),    ("internal energy", e, e_ex)]

for a, (name, num, ex) in zip(ax.ravel(), panels):
    a.plot(d["x"], ex,  "k-", lw=1.2, label="exact")
    a.plot(d["x"], num, "o", ms=3, mfc="none", color="tab:red", label="Godunov")
    a.set_title(name); a.grid(alpha=0.3)
    
ax[0, 0].legend()
fig.tight_layout()
fig.savefig("godunov.png", dpi=130)
plt.show()