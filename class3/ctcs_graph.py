#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# User config (match Fortran)
# -----------------------
fname = "advection_ctcs_RAW_filter.txt"

# spatial grid (match your Fortran)
x0 = 0.0
x1 = 500.0
dx = 0.1

# time sampling used in Fortran saving
tp = 200.0   # you save every ~tp (200)
t0 = 0.0     # initial time

# curves to plot (by snapshot index k)
snap_k = [0, 1, 2, 3, 4, 5] 

# -----------------------
# Load
# -----------------------
A = np.loadtxt(fname)  # shape: (nx, nsaved)
nx, nsaved = A.shape

x = x0 + dx * np.arange(nx)

# clip snapshot list to what's available
snap_k = [k for k in snap_k if 0 <= k < nsaved]
if not snap_k:
    snap_k = [0, nsaved // 2, nsaved - 1]

# physical times for each saved column
t_cols = t0 + tp * np.arange(nsaved)

# -----------------------
# Plot (style like your example)
# -----------------------
plt.figure(figsize=(10.5, 6.0))

for k in snap_k:
    t_phys = t_cols[k]
    plt.plot(x, A[:, k], linewidth=1.6, label=rf"$\phi^{{({int(t_phys)})}}$")

plt.title("Leapfrog scheme")

plt.xlabel("x", fontsize=12)
plt.ylabel(r"$\phi$", fontsize=12)

plt.grid(True, alpha=0.25)

# Make y-limits robust (ignore extreme outliers if any)
ymin, ymax = np.percentile(A[:, snap_k], [1, 99])
pad = 0.15 * (ymax - ymin + 1e-12)
plt.ylim(ymin - pad, ymax + pad)

plt.xlim(x.min(), x.max())

plt.legend(loc="upper right", framealpha=0.9, fontsize=11)
plt.tight_layout()
plt.show()

