import xarray as xr
import matplotlib.pyplot as plt

# file
ds = xr.open_dataset("advection_lagrangian.nc")

x = ds["x"]
time = ds["time"]

phi_lin = ds["phi_linear"]
phi_cub = ds["phi_cubic"]

# linear
plt.figure()

for t in time.values:
    phi_lin.sel(time=t).plot(x="x", label=f"t = {int(t)} s")

plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.title("Semi-Lagrangian Scheme (Linear Interpolation)")
plt.legend()
plt.tight_layout()

plt.show()
plt.close()

# cubic
plt.figure()

for t in time.values:
    phi_cub.sel(time=t).plot(x="x", label=f"t = {int(t)} s")

plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.title("Semi-Lagrangian Scheme (Cubic Interpolation)")
plt.legend()
plt.tight_layout()

plt.show()
plt.close()

