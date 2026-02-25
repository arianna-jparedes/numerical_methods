import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("euler_heun.txt")
x = data[:,0]
err_eu = data[:,1]
err_he = data[:,2]

plt.figure()
plt.plot(x, err_eu, label="Euler error")
plt.plot(x, err_he, label="Heun error")
plt.xlabel("x")
plt.ylabel("y - y_exact")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("errors.png", dpi=200)
plt.show()

print("max|error| Euler:", np.max(np.abs(err_eu)))
print("max|error| Heun :", np.max(np.abs(err_he)))

