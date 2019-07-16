import numpy as np
import matplotlib.pyplot as plt

# Code Units
A = 0.1
for N in np.power(8 * [2], range(9,17)):

    x = np.linspace(0, 1, N, endpoint=False)
    delta = -A * np.cos(2*np.pi * x)

    np.savetxt('./delta_%d.txt' % N, delta)
