import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def momentum(file, p=None):
    """
    Given an HDF5 file this function returns time-series of the total momentum

    Input: HDF5 file handle
    Returns: tuple (t, nump array (len(t))) where t is either code time if
    cosmo model is static or redshift value
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    entropies = file["Momentum"]
    tau, z  = [], []
    momentum_data = np.empty(len(entropies.values()))
    for i,entropy in enumerate(entropies.values()):
        tau.append(entropy.attrs["tau"][0])
        z.append(entropy.attrs["z"][0])
        momentum_data[i] = np.array(entropy)
    t = np.array(tau) if p["Cosmology"]["model"] == 1 else np.array(z)

    return (t, momentum_data)

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)

    t, P = momentum(file,p)
    fig,ax = plt.subplots()

    t_label = r"$\tau$" if p["Cosmology"]["model"] == 1 else r"$z$"


    ax.plot(t, P, label="P")
    ax.legend(*ax.get_legend_handles_labels())
    ax.set(xlabel=t_label)
    plt.show()
