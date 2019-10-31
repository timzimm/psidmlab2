import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def wavefunction(file, time, p=None):
    """
    Given an HDF5 file and a list of time points this function returns 
    a tuple of numpy arrays holding |psi|^2, Re(psi), Im(psi), Arg(psi)

    Input: HDF5 file handle, list of timepoints
    Returns: tuple of arrays (len(time) x N)
    """
    # Get simultion parameters and parse string as json
    if p is None:
        param_string = file['/'].attrs['parameters'][0].decode("ascii")
        p = json.loads(param_string)

    time_file=None
    psis = file["WaveFunction"]
    ds_names=list(psis.keys())

    if p["Cosmology"]["model"] == 1:
        time_file = np.array([psi.attrs['tau'][0] for psi in psis.values()])
    else:
        time_file = np.array([psi.attrs['z'][0] for psi in psis.values()])

    psis = [np.array(psis[ds_names[np.argmin(np.abs(time_file - t))]]) for t
                       in time]

    return psis

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    time =np.array(sys.argv[2:],dtype=float)

    # Get simultion parameters and parse string as json
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)
    psis = wavefunction(file, time, p)

    fig, axs = plt.subplots(4, 1, sharex=True)


    for psi, t in zip(psis, time):
        psi = np.ravel(psi)
        label = ""
        if(p["Cosmology"]["model"] == 1):
            label=r"$\tau = %.4f$" % t
        if(p["Cosmology"]["model"] == 0):
            label=r"$ z = %.4f$" % t
        axs[0].plot(np.abs(psi)**2, label=label)
        axs[1].plot(np.real(psi))
        axs[2].plot(np.imag(psi))
        axs[3].plot(np.angle(psi))

    axs[0].set(ylabel=r"$|\psi|^2$")
    axs[1].set(ylabel=r"$Re(\psi)$")
    axs[2].set(ylabel=r"$Im(\psi)$")
    axs[3].set(ylabel=r"$\Phi(\psi)$")

    for ax in axs:
        ax.grid()

    plt.legend(*axs[0].get_legend_handles_labels(), loc=0)
    plt.show()
