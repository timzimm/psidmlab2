import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import seaborn as sns
from matplotlib.colors import LogNorm 
import sys

file = h5py.File(sys.argv[1], 'r')
time = np.float64(sys.argv[2])
print(time)

param_string = file['/'].attrs['parameters'][0].decode("ascii")
p = json.loads(param_string)

L = int(p["Simulation"]["L"])
N = int(p["Simulation"]["N"])
sigma_x = float(p["Observables"]["PhaseSpaceDistribution"]["sigma_x"])

phasedists = file["PhaseSpaceDistribution"]

fig,ax = plt.subplots()

dist = None

if p["Cosmology"]["model"] == 0:
    dist = [phasedists[i] for i in phasedists if
                np.isclose(phasedists[i].attrs['tau'][0], time)][0]
else:
    dist = [phasedists[i] for i in phasedists if
                np.isclose(phasedists[i].attrs['z'][0], time)][0]


cm = sns.diverging_palette(240, 10, n=9, as_cmap=True)
lvl = np.logspace(-2,0,6)

tau = dist.attrs["tau"][0]

if bool(dist.attrs["transpose"][0]):
    dist = np.transpose(np.array(dist))

dist = np.array(dist)[46*N//100:54*N//100]
dist /= np.max(dist)

x = np.linspace(0, L, N)
u = np.linspace(-2*np.pi/L, 2*np.pi/L, N)[46*N//100:54*N//100]
X, U = np.meshgrid(x,u)

print("Creating pcolormesh...")
cbar = ax.contourf(X, U, dist, levels=lvl,norm=LogNorm(),cmap=cm)
# cbar = ax.contourf(X, U, np.clip(dist,0,0.001))
fig.colorbar(cbar, ax=ax)

plt.title(r"Phasespace Distribution ($\sigma_x = %.3f$)" % sigma_x)
plt.show()
