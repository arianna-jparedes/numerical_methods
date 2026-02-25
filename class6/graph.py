import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset("diffusion_tdma.nc")

x = ds["x"]
time = ds["time"]

phi = ds["phi"].transpose("time","x")

for t in time.values:
    phi.sel(time=t).plot(x="x", label=f"t = {int(t)} s")

plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.legend()
plt.title("Diffusion Metal Rod (TDMA method)")
plt.tight_layout()
plt.show()

