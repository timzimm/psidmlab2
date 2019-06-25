import numpy as np
import matplotlib.pyplot as plt

# Code Units
L = 20 
omega = 1;
N = 1024
cols = 2

pot = np.empty(N*cols).reshape(N, cols)
pot[:,0] = np.linspace(0, L, N, endpoint=False)
pot[:,1] = 0.5 * omega ** 2 *(pot[:,0] - L/2)**2

psi = np.empty(N*cols).reshape(N, cols)
psi[:,0] = (omega/np.pi)**0.25 * np.sqrt(L)* np.exp(-0.5 * omega * (pot[:,0] - L/4)**2)
psi[:,1] = 0

# Second Wavefunction
psi[:,0] += np.flip(psi[:,0]) 

#Check normalization
print(L/N * np.sum(psi[:,0]**2))

plt.plot(pot[:,0], pot[:,1])
plt.plot(pot[:,0], psi[:,0])
plt.show()

np.savetxt('./potential.txt', pot, header="x \t U")
np.savetxt('./psi.txt', psi, header="Re(psi) \t Im(psi)")
