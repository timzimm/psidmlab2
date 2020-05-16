import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import seaborn as sns
from matplotlib.colors import LogNorm 
import sys

def phasespace(file, time, p=None):
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    time_file=None
    fs = file["PhaseSpaceDistribution"]
    ds_names=list(fs.keys())

    if p["Cosmology"]["model"] == 1:
        time_file = np.array([f.attrs['tau'][0] for f in fs.values()])
    else:
        time_file = np.array([f.attrs['z'][0] for f in fs.values()])

    fs = [np.array(fs[ds_names[np.argmin(np.abs(time_file - t))]]) for t
                       in time]

    return fs

if __name__ == "__main__":
    p = None
    dists = None
    time = np.float64(sys.argv[2:])
    with h5py.File(sys.argv[1], 'r') as file:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)
        dists = np.array(phasespace(file, time, p))

    L = int(p["Simulation"]["L"])
    N = int(p["Simulation"]["N"])
    sigma_x = float(p["Observables"]["PhaseSpaceDistribution"]["sigma_x"])

    fig,ax = plt.subplots()
    cm = sns.diverging_palette(240, 10, n=9, as_cmap=True)

    for i in range(len(time)):
        dist = dists[i,:,:]
        dist /= np.max(dist)
        print(np.all(dist>0))

        patch = np.empty((2,2))
        patch[0,0] = -L/2
        patch[0,1] = -patch[0,0] - L/N
        patch[1,0] = -N*np.pi/L
        patch[1,1] = -patch[1,0] - 2*np.pi/L
        # patch = np.array(p["Observables"]["PhaseSpaceDistribution"]["patch"])
        # Lx = patch[0,1] - patch[0,0]
        # Lk = patch[1,1] - patch[1,0]

        x = np.linspace(patch[0,0], patch[0,1], dist.shape[1])
        k = np.linspace(patch[1,0], patch[1,1], dist.shape[0])

        X, K = np.meshgrid(x,k)
        dist = np.clip(dist,0.001,1)
        dist[dist <= 0.001] = None
        print(dist.shape)
        cbar = ax.contourf(X, K, dist, 
                           levels=np.logspace(-3,0,100), norm=LogNorm(), 
                           cmap="afmhot")
        fig.colorbar(cbar, ax=ax)

        plt.title(r"Phasespace Distribution ($\sigma_x = %.3f$)" % sigma_x)
        plt.show()
