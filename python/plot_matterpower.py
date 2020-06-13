import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys
def matter_power_spectrum(file, time, p=None):
    """
    Given an HDF5 file and a list of time points this function returns 
    a numpy array holding the matter power spectrum at the requested time instances
    (or as close to them as possible)

    Input: HDF5 file handle, list of timepoints
    Returns: numpy array len(time) x N
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    time_file=None
    Ps = file["MatterPowerSpectrum"]
    ds_names=list(Ps.keys())

    if p["Cosmology"]["model"] == 1:
        time_file = np.array([P.attrs['tau'][0] for P in Ps.values()])
    else:
        time_file = np.array([P.attrs['z'][0] for P in Ps.values()])

    Ps = [np.array(Ps[ds_names[np.argmin(np.abs(time_file - t))]]) for t
                       in time]
    return np.array(Ps)

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    time = np.array(sys.argv[2:],dtype=float)
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)
    Ps = matter_power_spectrum(file, time, p)

    fig,ax = plt.subplots()

    for i,t in enumerate(time):
        if(p["Cosmology"]["model"] == 0):
            ax.loglog(Ps[i,:], label=r"$\tau = %.2f$" % t)
        if(p["Cosmology"]["model"] == 1):
            ax.loglog(Ps[i,:], label=r"$ z = %.2f$" % t)

    ax.set(xlim=[1, int(p["Domain"]["N"])//2 + 1], ylabel=r"$P$")
    ax.legend(*ax.get_legend_handles_labels())
    ax.grid()
    plt.show()
