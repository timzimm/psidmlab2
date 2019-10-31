import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def energy(file, p=None):
    """
    Given an HDF5 file this function returns time-series of kinnetic, potential
    ans total energy, as well as the discrepancy of the virial theorem

    Input: HDF5 file handle
    Returns: tuple (T, V, <T> - <virial>, t) where t is either code time if
    cosmo model is static or redshift value
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    energies = file["Energy"]
    tau, z  = [], []
    energy_data = np.empty((len(energies.values()), 4))
    for i,energy in enumerate(energies.values()):
        tau.append(energy.attrs["tau"][0])
        z.append(energy.attrs["z"][0])
        energy_data[i,:] = np.array(energy)

    t = np.array(tau) if p["Cosmology"]["model"] == 1 else np.array(z)
    T = energy_data[:,0]
    V = energy_data[:,1]
    Etot = energy_data[:,2]
    virial = energy_data[:,3]

    #Time average with adaptive dt
    T_mean = np.divide(np.cumsum(np.multiply(np.diff(tau), T[1:])),
                       tau[1:])
    virial_mean = np.divide(np.cumsum(np.multiply(np.diff(tau), virial[1:])),
                            tau[1:])

    return T, V, Etot, T_mean, virial_mean, t

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)

    T, V, Etot, T_mean, virial_mean, t = energy(file, p)
    delta_virial = np.abs(2*T_mean - virial_mean)
    fig,ax = plt.subplots(nrows=2, ncols=1, sharex=True)

    t_label = r"$\tau$" if p["Cosmology"]["model"] == 0 else r"$a$"


    ax[0].plot(t, T, label="T")
    ax[0].plot(t, V, label="V")
    ax[0].plot(t, Etot, label="Etot")
    ax[1].loglog(t[1:], delta_virial, 
            label=r"$2\langle T \rangle - \langle x \partial_x V\rangle$")

    for axis in ax:
        axis.set(ylabel=r"$E$")
        axis.grid()
    ax[0].legend(*ax[0].get_legend_handles_labels())
    ax[1].legend(*ax[1].get_legend_handles_labels())
    ax[1].set(xlabel=t_label)
    plt.show()
