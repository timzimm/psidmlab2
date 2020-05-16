import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def entropy_density(file, time, p=None):
    """
    Given an HDF5 file and a list of time points this function returns 
    a numpy array holding the entropy density at the requested time instances
    (or as close to them as possible)

    Input: HDF5 file handle, list of timepoints
    Returns: numpy array len(time) x N
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    time_file=None
    sdenisty = file["EntropyDensity"]
    ds_names=list(sdenisty.keys())

    if p["Cosmology"]["model"] == 1:
        time_file = np.array([s.attrs['tau'][0] for s in sdenisty.values()])
    else:
        time_file = np.array([s.attrs['z'][0] for s in sdenisty.values()])

    sdenisty = [np.array(sdenisty[ds_names[np.argmin(np.abs(time_file - t))]]) for t
                       in time]
    return sdenisty

if __name__ == "__main__":
    time = np.array(sys.argv[2:],dtype=float)
    p = None
    deltas = None
    with h5py.File(sys.argv[1], 'r') as file:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)
        sdensity = entropy_density(file, time, p)

    fig,ax = plt.subplots()

    for i,t in enumerate(time):
        if(p["Cosmology"]["model"] == 0):
            ax.plot(sdensity[i,:], label=r"$\tau = %.2f$" % t)
        if(p["Cosmology"]["model"] == 1):
            ax.plot(sdensity[i,:], label=r"$ z = %.2f$" % t)

    ax.set(xlim=[0, p["Simulation"]["N"]], ylabel=r"$s$")
    ax.legend(*ax.get_legend_handles_labels())
    ax.grid()
    plt.show()
