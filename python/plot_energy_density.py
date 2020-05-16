import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def energy_densities(file, time, p=None):
    """
    Given an HDF5 file and a list of time points this function returns 
    a list of numpy arrays for each requested time point.

    Input: HDF5 file handle, list of timepoints
    Returns: list (len(time)) of arrays (3xN) with
    array[0,:] = kinetic energy density
    array[1,:] = potential energy density
    array[2,:] = virial density
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    time_file=None
    Es = file["EnergyDensity"]
    ds_names=list(Es.keys())

    if p["Cosmology"]["model"] == 1:
        time_file = np.array([E.attrs['tau'][0] for E in Es.values()])
    else:
        time_file = np.array([E.attrs['z'][0] for E in Es.values()])

    Elist = []
    for t in time:
        Elist.append(np.array(Es[ds_names[np.argmin(np.abs(time_file - t))]]))

    return Elist

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    time =np.array(sys.argv[2:],dtype=float)

    # Get simultion parameters and parse string as json
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)
    Es = energy_densities(file, time, p)
    N = p["Simulation"]["N"]
    L = p["Simulation"]["L"]
    x = np.linspace(-L/2,L/2,N, endpoint=False)

    fig, axs = plt.subplots(3, 1, sharex=True)

    for E, t in zip(Es, time):
        label = ""
        if(p["Cosmology"]["model"] == 1):
            label=r"$\tau = %.4f$" % t
        if(p["Cosmology"]["model"] == 0):
            label=r"$ z = %.4f$" % t
        axs[0].plot(x,E[0,:], label=label)
        axs[1].plot(x,E[1,:])
        axs[2].plot(x,E[2,:])

    axs[0].set(ylabel=r"$T$")
    axs[1].set(ylabel=r"$V$")
    axs[2].set(ylabel=r"$Virial$")

    for ax in axs:
        ax.grid()

    plt.legend(*axs[0].get_legend_handles_labels(), loc=0)
    plt.show()
