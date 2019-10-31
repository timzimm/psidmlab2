import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def flux(file, time, p=None):
    """
    Given an HDF5 file and a list of time points this function returns 
    a numpy array holding the particle flux

    Input: HDF5 file handle, list of timepoints
    Returns: numpy array (len(time) x N)
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    time_file=None
    fluxes = file["ParticleFlux"]
    ds_names=list(fluxes.keys())

    if p["Cosmology"]["model"] == 1:
        time_file = np.array([flux.attrs['tau'][0] for flux in fluxes.values()])
    else:
        time_file = np.array([flux.attrs['z'][0] for flux in fluxes.values()])

    fluxes = [np.array(fluxes[ds_names[np.argmin(np.abs(time_file - t))]]) for t
                       in time]
    return np.array(fluxes)

if __name__ == "__main__":
    time =np.array(sys.argv[2:],dtype=float)
    file = h5py.File(sys.argv[1], 'r')
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)
    fluxes = flux(file, time, p)

    fig,ax = plt.subplots()

    for i, t in zip(range(psi2.shape[0]), time):
        if(p["Cosmology"]["model"] == 0):
            ax.plot(fluxes[i,:], label=r"$\tau = %.2f$" % t)
        if(p["Cosmology"]["model"] == 1):
            ax.plot(fluxes[i,:], label=r"$ z = %.2f$" % t)

    ax.set(xlim=[0, p["Simulation"]["N"]], ylabel=r"$j$")
    ax.legend(*ax.get_legend_handles_labels())
    ax.grid()
    plt.show()
