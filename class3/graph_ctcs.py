import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset("advection_ctcs.nc")

x = ds["x"]
time = ds["time"]
phi = ds["phi nonfilter"]
phi_ra = ds["phi RA filter"]
phi_raw = ds["phi RAW filter"]

# Plot snapshots

for t in time.values:
    phi.sel(time=t).plot(x="x", label=f"t = {int(t)} s")

plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.legend()
plt.title("Linear advection (non filter)")
plt.tight_layout()
plt.savefig("./adv_nonfilter.png", dpi=200)
plt.show()

for t in time.values:
    phi_ra.sel(time=t).plot(x="x", label=f"t = {int(t)} s")

plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.legend()
plt.title("Linear advection (RA filter)")
plt.tight_layout()
plt.savefig("./adv_rafilter.png", dpi=200)
plt.show()

for t in time.values:
    phi_raw.sel(time=t).plot(x="x", label=f"t = {int(t)} s")

plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.legend()
plt.title("Linear advection (RAW filter)")
plt.tight_layout()
plt.savefig("./adv_rawfilter.png", dpi=200)
plt.show()



