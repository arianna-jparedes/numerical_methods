import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# file
ds = xr.open_dataset("gravity_wave.nc")
t_req = np.arange(0, 2001, 200)
sel = ds.sel(time=t_req, method="nearest")

# plot of psi
fig, ax = plt.subplots()
for t in sel.time.values:
    ax.plot(sel["x"].values, sel["p"].sel(time=t).values, label=f"t={float(t):.0f}s")
ax.set_xlabel("x (m)")
ax.set_ylabel("p (Φ)  [m² s⁻²]")
ax.set_title("Gravity wave: Φ(x,t) Method 1")
ax.legend(ncol=1, fontsize=8)
plt.tight_layout()
plt.show()

# plot of u 
fig, ax = plt.subplots()
for t in sel.time.values:
    ax.plot(sel["x"].values, sel["u"].sel(time=t).values, label=f"t={float(t):.0f}s")
ax.set_xlabel("x (m)")
ax.set_ylabel("u  [m s⁻¹]")
ax.set_title("Gravity wave: u(x,t) Method 1")
ax.legend(ncol=1, fontsize=8)
plt.tight_layout()
plt.show()

# plot of psi
fig, ax = plt.subplots()
for t in sel.time.values:
    ax.plot(sel["x"].values, sel["p2"].sel(time=t).values, label=f"t={float(t):.0f}s")
ax.set_xlabel("x (m)")
ax.set_ylabel("p (Φ)  [m² s⁻²]")
ax.set_title("Gravity wave: Φ(x,t) Method 2")
ax.legend(ncol=1, fontsize=8)
plt.tight_layout()
plt.show()

# plot of u 
fig, ax = plt.subplots()
for t in sel.time.values:
    ax.plot(sel["x"].values, sel["u2"].sel(time=t).values, label=f"t={float(t):.0f}s")
ax.set_xlabel("x (m)")
ax.set_ylabel("u  [m s⁻¹]")
ax.set_title("Gravity wave: u(x,t) Method 2")
ax.legend(ncol=1, fontsize=8)
plt.tight_layout()
plt.show()


