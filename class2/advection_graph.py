import numpy as np
import matplotlib.pyplot as plt

dx = 0.1
x0 = 0.0

data = np.loadtxt("advection.txt")

nx, nsnap = data.shape
x = x0 + dx * np.arange(nx)

tp = 200
times = np.arange(nsnap) * tp

for k in range(nsnap):
    plt.plot(x, data[:,k], label=f"t = {times[k]} s")

plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.legend()
plt.show()

