import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def scalefactor(file, p=None):
    """
    Given an HDF5 file returns the numerical values of code time and the scalefactor at the
    times observables were computed
    """
    if p is None:
        # Get simultion parameters and parse string as json
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    first_observable = file[list(file.keys())[0]]
    tau = np.empty(len(first_observable))
    a = np.empty(len(first_observable))
    for i,obs in enumerate(first_observable.values()):
        tau[i] = obs.attrs["tau"][0]
        a[i] = obs.attrs["a"][0]
    return tau,a

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)
    tau, a = scalefactor(file, p)

    fig,ax = plt.subplots()
    ax.plot(tau, a)
    ax.set(ylabel=r"$a$", xlabel=r"$\tau$")
    ax.grid()
    plt.show()
