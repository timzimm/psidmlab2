import numpy as np
import matplotlib.pyplot as plt

# Code Units
L = 1
beta = 128/L
N = 16384

for N in np.power(8 * [2], range(9,17)):
    x = np.linspace(-0.5, 0.5, N, endpoint=False)
    delta = np.tanh(beta*(x+L/20)) + np.tanh(beta*(L/20-x))
    delta -= L/N * np.sum(delta)
    plt.plot(x, delta)
    plt.show()
    np.savetxt('./delta_%d.txt' % N, delta)
