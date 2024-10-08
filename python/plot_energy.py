import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def energy(file, p=None):
    """
    Given an HDF5 file this function returns time-series of kinetic, potential
    and the virial.

    Input: HDF5 file handle
    Returns: tuple (t, nump array (len(t) x 3)) where t is either code time if
    cosmo model is static or redshift value
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    energies = file["Energy"]
    tau, z  = [], []
    energy_data = np.empty((len(energies.values()), 3))
    for i,energy in enumerate(energies.values()):
        tau.append(energy.attrs["tau"][0])
        z.append(energy.attrs["z"][0])
        energy_data[i,:] = np.array(energy)

    t = np.array(tau) if p["Cosmology"]["model"] == 1 else np.array(z)

    return (t, energy_data)

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)

    t, Es = energy(file,p)
    fig,ax = plt.subplots()

    t_label = r"$\tau$" if p["Cosmology"]["model"] == 1 else r"$z$"


    ax.plot(t, Es[:,0], label="T")
    ax.plot(t, Es[:,1], label="V")
    ax.plot(t, Es[:,0]+Es[:,1], label=r"$T+V$")
    ax.legend(*ax.get_legend_handles_labels())
    ax.set(xlabel=t_label)
    plt.show()
