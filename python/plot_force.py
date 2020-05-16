import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def force(file, time, p=None):
    """
    Given an HDF5 file and a list of time points this function returns 
    a numpy array holding the forces at the requested time instances
    (or as close to them as possible)

    Input: HDF5 file handle, list of timepoints
    Returns: numpy array len(time) x 2 x N
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    time_file=None
    forces = file["Force"]
    ds_names=list(forces.keys())

    if p["Cosmology"]["model"] == 1:
        time_file = np.array([f.attrs['tau'][0] for f in forces.values()])
    else:
        time_file = np.array([f.attrs['z'][0] for f in forces.values()])

    forces = [np.array(forces[ds_names[np.argmin(np.abs(time_file - t))]]) for t
                       in time]
    return forces

if __name__ == "__main__":
    time = np.array(sys.argv[2:],dtype=float)
    p = None
    deltas = None
    with h5py.File(sys.argv[1], 'r') as file:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)
        forces = force(file, time, p)

    fig,ax = plt.subplots(nrows=2)

    for i,t in enumerate(time):
        if(p["Cosmology"]["model"] == 0):
            ax[0].plot(forces[i,0,:], label=r"$\tau = %.2f$" % t)
            ax[1].plot(forces[i,1,:], label=r"$\tau = %.2f$" % t)
        if(p["Cosmology"]["model"] == 1):
            ax[0].plot(forces[i,0,:], label=r"$ z = %.2f$" % t)
            ax[1].plot(forces[i,1,:], label=r"$ z = %.2f$" % t)

    ax.set(xlim=[0, p["Domain"]["N"]], ylabel=r"$F$")
    ax.legend(*ax.get_legend_handles_labels())
    ax.grid()
    plt.show()
