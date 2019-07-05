import numpy as np
import matplotlib.pyplot as plt

# Code Units
L = 1 
A = 0.1
N = 2**10

x = np.linspace(0, L, N, endpoint=False)
dx = x[1] - x[0]

delta = -A * np.cos(2*np.pi/L * x)

print(np.sum(delta))

plt.plot(x, delta)
plt.show()

np.savetxt('./delta.txt', delta)
