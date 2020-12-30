import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import json
import sys
from copy import deepcopy
from pathlib import Path
import argparse


def a_of(z):
    return 1/(1+z)


def z_of(a):
    return 1/a - 1


def dtda(t, a, om, ol):
    om_at_a = om/(om + ol * a*a*a)
    return np.sqrt(1.5 * om_at_a / (a*a*a))


parser = argparse.ArgumentParser(
    description="Chunk psidm json into multiple parts")
parser.add_argument("chunks",
                    help="number of equal size chunks or list of chunk start times"
                    )
parser.add_argument("file", help="json file to process")

args = parser.parse_args()

with open(args.file, "r") as file:
    config = json.load(file)
    cosmo_model = not bool(config["Cosmology"]["model"])
    om = float(config["Cosmology"]["omega_m0"]) if cosmo_model else None

    # Checkpoints might be code time t or redshift z, depending on the value of
    # cosmo_model
    checkpoints = set()
    for observables in config["Observables"].items():
        compute_at = observables[1]["compute_at"]
        if compute_at != [-1]:
            checkpoints.update(compute_at)
    checkpoints = sorted(list(checkpoints), reverse=cosmo_model)


def t_of_checkpoint(checkpoint): return checkpoint


def checkpoint_of_t(t): return t


if cosmo_model:
    # Construct function t(z)
    amin = a_of(checkpoints[0])
    amax = a_of(checkpoints[-1])
    a_grid = np.linspace(amin, amax, 1000)
    t = np.asarray(odeint(dtda, 0, a_grid, args=(om, 1-om))).flatten()
    t_of_checkpoint = interp1d(z_of(a_grid), t)
    checkpoint_of_t = interp1d(t, z_of(a_grid))

args.chunks = args.chunks.split(",")
last_checkpoint = []
if len(args.chunks) == 1:
    chunk_count = int(args.chunks[0])

    # Assure equal runtime by dividing total code time
    delta_t = np.abs(t_of_checkpoint(checkpoints[0]) -
                     t_of_checkpoint(checkpoints[-1]))/chunk_count

    # Adjust accuracy of chunk sizes to accuracy of exsisting checkpoints.
    # Not really nececssary, just cosmetics.
    decimals = max(0, np.amax([f"{c}"[::-1].find('.') for c in checkpoints]))
    last_checkpoint = [checkpoint_of_t([(i+1) * delta_t])[0]
                       for i in range(chunk_count)]
    last_checkpoint = np.round(last_checkpoint, decimals=decimals).tolist()
else:
    args.chunks = args.chunks[:-1] if args.chunks[-1] == "" else args.chunks
    last_checkpoint = [int(cp) for cp in args.chunks]

last_t = t_of_checkpoint(last_checkpoint)

# Add sentinal checkpoints to wave function to guarantee equal chunk runtimes
psi_checkpoints = config["Observables"]["WaveFunction"]["compute_at"]
psi_checkpoints += last_checkpoint

if psi_checkpoints[0] < 0:
    psi_checkpoints = psi_checkpoints[1:] + checkpoints

# Remove entries occcuring more than once
psi_checkpoints = list(set(psi_checkpoints))
psi_checkpoints = sorted(psi_checkpoints, reverse=cosmo_model)
config["Observables"]["WaveFunction"]["compute_at"] = psi_checkpoints

output_path = Path(config["General"]["output_file"])
json_path = Path(args.file)
json_dir = json_path.parent

for i, next_t in enumerate(last_t):
    config_chunk = deepcopy(config)
    config_obs = config_chunk["Observables"]
    config_gen = config_chunk["General"]
    config_ic = config_chunk["Initial Conditions"]

    # Store in same HDF5 file as previous chunk
    config_gen["output_file"] = str(output_path)

    # Remove checkpoints in each observable excceeding next_t
    for observable in config_obs.items():
        compute_at = np.array(sorted(observable[1]["compute_at"],
                                     reverse=cosmo_model))
        if len(compute_at) != 0 and compute_at[0] >= 0:
            observable[1]["compute_at"] = \
                compute_at[t_of_checkpoint(compute_at) <= next_t].tolist()

    if i > 0:
        # use last output file as initial condition
        config_ic["ic_type"] = 3
        config_ic["source_file"] = str(output_path)

    json_new = json_path.stem + f"_chunk_{i}" + json_path.suffix
    with open(json_dir/json_new, 'w') as file:
        json.dump(config_chunk, file, indent=4)
