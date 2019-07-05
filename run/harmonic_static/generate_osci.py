import numpy as np
import matplotlib.pyplot as plt

# Code Units
L = 20 
omega = 1;
N = 1024

x = np.linspace(0, L, N, endpoint=False)
pot = 0.5 * omega ** 2 *(x - L/2)**2

# First Wavefunction
psi = np.exp(-0.5 * omega * (x - L/4)**2)
# Second Wavefunction
psi += np.flip(psi)
# Normalization
psi /= np.sqrt(L/N * np.sum(psi**2))
print(L/N * np.sum(psi**2))

delta = L * psi**2 - 1

#Check normalization
print(np.sum(delta))

plt.plot(x, pot)
plt.plot(x, delta)
plt.show()

np.savetxt('./potential.txt', pot)
np.savetxt('./delta.txt', delta)
