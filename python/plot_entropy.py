import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def entropy(file, p=None):
    """
    Given an HDF5 file this function returns time-series of the entropy

    Input: HDF5 file handle
    Returns: tuple (t, nump array (len(t))) where t is either code time if
    cosmo model is static or redshift value
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    entropies = file["Entropy"]
    tau, z  = [], []
    entropy_data = np.empty(len(entropies.values()))
    for i,entropy in enumerate(entropies.values()):
        tau.append(entropy.attrs["tau"][0])
        z.append(entropy.attrs["z"][0])
        entropy_data[i] = np.array(entropy)
    t = np.array(tau) if p["Cosmology"]["model"] == 1 else np.array(z)

    #Time average with adaptive dt
    # T_mean = np.divide(np.cumsum(np.multiply(np.diff(tau), T[1:])),
    #                    tau[1:])
    # virial_mean = np.divide(np.cumsum(np.multiply(np.diff(tau), virial[1:])),
    #                         tau[1:])

    # return T, V, Etot, T_mean, virial_mean, t
    return (t, entropy_data)

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)

    t, Es = energy(file,p)
    # T, V, Etot, T_mean, virial_mean, t = energy(file, p)
    # delta_virial = np.abs(2*T_mean - virial_mean)
    fig,ax = plt.subplots()

    t_label = r"$\tau$" if p["Cosmology"]["model"] == 0 else r"$z$"


    ax.plot(t, Es[:,0], label="T")
    ax.plot(t, Es[:,1], label="V")
    ax.plot(t, Es[:,0]+Es[:,1], label=r"$T+V$")
    # ax.plot(t, Es[:,2], label="virial")
    ax.legend(*ax.get_legend_handles_labels())
    ax.set(xlabel=t_label)
    # ax[1].loglog(t[1:], delta_virial, 
    #         label=r"$2\langle T \rangle - \langle x \partial_x V\rangle$")

    # for axis in ax:
    #     axis.set(ylabel=r"$E$")
    #     axis.grid()
    # ax[0].legend(*ax[0].get_legend_handles_labels())
    # ax[1].legend(*ax[1].get_legend_handles_labels())
    # ax[1].set(xlabel=t_label)
    plt.show()
