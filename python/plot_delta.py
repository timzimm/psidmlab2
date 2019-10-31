import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def density_contrast(file, time, p=None):
    """
    Given an HDF5 file and a list of time points this function returns 
    a numpy array holding the density contrast at the requested time instances
    (or as close to them as possible)

    Input: HDF5 file handle, list of timepoints
    Returns: numpy array len(time) x N
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    time_file=None
    deltas = file["DensityContrast"]
    ds_names=list(deltas.keys())

    if p["Cosmology"]["model"] == 1:
        time_file = np.array([delta.attrs['tau'][0] for delta in deltas.values()])
    else:
        time_file = np.array([delta.attrs['z'][0] for delta in deltas.values()])

    deltas = [np.array(deltas[ds_names[np.argmin(np.abs(time_file - t))]]) for t
                       in time]
    return np.array(deltas)

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    time = np.array(sys.argv[2:],dtype=float)
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)
    deltas = density_contrast(file, time, p)

    fig,ax = plt.subplots()

    for i,t in enumerate(time):
        if(p["Cosmology"]["model"] == 0):
            ax.plot(deltas[i,:], label=r"$\tau = %.2f$" % t)
        if(p["Cosmology"]["model"] == 1):
            ax.plot(deltas[i,:], label=r"$ z = %.2f$" % t)

    ax.set(xlim=[0, p["Simulation"]["N"]], ylabel=r"$\delta$")
    ax.legend(*ax.get_legend_handles_labels())
    ax.grid()
    plt.show()
