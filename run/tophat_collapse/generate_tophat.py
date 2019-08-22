import numpy as np
import matplotlib.pyplot as plt

# Code Units
L = 1
beta = 10/L
Ns = [2**14, 2**15]

for N in Ns:
    for width in [0.01, 0.02, 0.04, 0.06, 0.08, 0.1]:
        x = np.linspace(-0.5, 0.5, N, endpoint=False)
        delta = np.tanh(beta/width*(x+0.5*width*L)) + np.tanh(beta/width*(0.5*width*L-x))
        delta -= L/N * np.sum(delta)
        plt.plot(x, delta)
        plt.show()
        np.savetxt('./ic/delta_%dkpc_%d.txt' % (width * 1000, N), delta)
