import xarray as xr
import matplotlib.pyplot as plt

# file
ds = xr.open_dataset("advection_diffusion.nc")

x = ds["x"]
time = ds["time"]

phi_lin = ds["phi"]

# linear
plt.figure()

for t in time.values:
    phi_lin.sel(time=t).plot(x="x", label=f"t = {int(t)} s")

plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.title("CTCS Scheme with diffusion")
plt.legend()
plt.tight_layout()

plt.show()
