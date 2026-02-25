import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset("advection.nc")

x = ds["x"]
time = ds["time"]
phi = ds["phi"]

# Plot snapshots
for t in time.values:
    phi.sel(time=t).plot(x="x", label=f"t = {int(t)} s")

plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.legend()
plt.tight_layout()
plt.savefig("./advection_nc.png", dpi=200)
plt.show()

