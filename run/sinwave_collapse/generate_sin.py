import numpy as np
import matplotlib.pyplot as plt

# Code Units
L = 1 
A = 0.1;
N = 1024
cols = 2

x = np.linspace(0, L, N, endpoint=False)

psi = np.empty(N*cols).reshape(N, cols)
psi[:,0] = np.sqrt(1.0/L * (A * np.sin(2*np.pi/L * x) + 1))
psi[:,1] = 0

plt.plot(x, psi[:,0])
plt.show()

np.savetxt('../run/psi_sin.txt', psi, header="Re(psi) \t Im(psi)")
