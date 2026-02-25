import xarray as xr
import matplotlib.pyplot as plt

# Read file
ncfile = "euler_heun.nc"
ds = xr.open_dataset(ncfile)

# Variables
x = ds["x"]
e_euler = ds["euler_error"]
e_heun  = ds["heun_error"]

# Graph
plt.figure(figsize=(10,6), dpi=150)
plt.plot(x, e_euler, linewidth=1.8, label="Euler scheme")
plt.plot(x, e_heun,  linewidth=1.8, label="Heun scheme")
plt.title(r"Error $(y_j - y(x_j))$ for ODE solver with $dx = 0.1$", fontsize=18)
plt.grid(True, alpha=0.25)
plt.xlim(0,10)
plt.ylim(-0.02,0.30)
plt.legend(loc="upper right", frameon=True)
plt.tight_layout()
plt.show()
