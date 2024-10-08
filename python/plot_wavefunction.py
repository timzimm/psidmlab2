import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import sys

def wavefunction(file, time=None, p=None):
    """
    Given an HDF5 file and a list of time points this function returns 
    a list of complex numpy arrays psi for each requested time point

    Input: HDF5 file handle, list of timepoints
    Returns: list (len(time)) of arrays (N)
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


    psis_complex = []
    if(time is not None):
        for t in time:
            psi_flat = np.ravel(psis[ds_names[np.argmin(np.abs(time_file - t))]])
            psis_complex.append(psi_flat[::2] + 1.0j * psi_flat[1::2])
    else:
        for name in ds_names:
            psi_flat = np.ravel(psis[name])
            psis_complex.append(psi_flat[::2] + 1.0j * psi_flat[1::2])


    return psis_complex

if __name__ == "__main__":
    file = h5py.File(sys.argv[1], 'r')
    time =np.array(sys.argv[2:],dtype=float)

    # Get simultion parameters and parse string as json
    param_string = file['/'].attrs['parameters'][0].decode("ascii")
    p = json.loads(param_string)
    psis = wavefunction(file, time, p)
    N = p["Domain"]["N"]
    L = p["Domain"]["L"]
    x = np.linspace(-L/2,L/2,N, endpoint=False)

    fig, axs = plt.subplots(5, 1, sharex=True)

    for psi, t in zip(psis, time):
        label = ""
        if(p["Cosmology"]["model"] == 1):
            label=r"$\tau = %.4f$" % t
        if(p["Cosmology"]["model"] == 0):
            label=r"$ z = %.4f$" % t
        axs[0].plot(x,np.abs(psi)**2, label=label)
        axs[1].plot(x,np.real(psi))
        axs[2].plot(x,np.imag(psi))
        phase = np.unwrap(np.angle(psi))
        axs[3].plot(x,phase)
        axs[4].plot(x,np.gradient(phase, L/N))

    axs[0].set(ylabel=r"$|\psi|^2$")
    axs[1].set(ylabel=r"$Re(\psi)$")
    axs[2].set(ylabel=r"$Im(\psi)$")
    axs[3].set(ylabel=r"$\Phi(\psi)$")
    axs[4].set(ylabel=r"$\partial_x\Phi(\psi)$")
    # axs[0].set_xlim(x[0],x[-1])

    for ax in axs:
        ax.grid()

    plt.legend(*axs[0].get_legend_handles_labels(), loc=0)
    plt.show()
