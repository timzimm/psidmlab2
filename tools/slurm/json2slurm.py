import numpy as np
import json
import time
import sys
import re
from pathlib import Path
from textwrap import dedent
from shutil import copyfile


def sec_to_hhmmss(time):
    hour, rem = divmod(time, 60 * 60)
    mins, rem = divmod(rem, 60)
    secs = rem % 60
    return "{0:02d}:{1:02d}:{2:02d}".format(hour, mins, secs)


def effify(non_f_str: str):
    """Converts string into f-string"""
    return eval(f'f"""{non_f_str}"""')


def sorted_nicely(l):
    """Sorts the given iterable in the way that is expected.

    Required arguments:
    l -- The iterable to be sorted.

    """
    def convert(text): return int(text) if text.isdigit() else text
    def alphanum_key(key): return [convert(c)
                                   for c in re.split("([0-9]+)", key)]
    return sorted(l, key=alphanum_key)


def memory(config):
    """Computes memory consumption of the simulation scenario specified in config"""
    # TODO: Make measurements, construct a mathematical model
    return 2000


def walltime(config):
    """Computes wall clock time of the simulation scenario specified in config"""
    # TODO: Make measurements, construct a mathematical model
    return 60 * 60 * 72


def topology(config):
    """Determines process execution topology i.e.:
    nodeN   -   number of computation nodes (required for distributed memory
                parallelization. Currently only nodeN = 1 supported
    taskN   -   number of tasks/processes available to job (i.e. MPI processes)
                Currently only taskN = 1 supported
    coresN  -   Number of cores per task. Used for shared memory parallelization
                Note:   bwunicluster2 allows hyperthreading. Thus if job is
                        exepcted to use N threads set coresN = 2*N and
                        OMP_NUM_THREADS=coresN/2
    """
    nodeN = 1
    taskN = 1
    coresN = 8
    return (nodeN, taskN, coresN)


def partition(config):
    """Determines best fitting queue of the simulation scenario specified in config"""
    partition_name = ""
    mem = memory(config)
    if mem > 180000:
        print(f"{mem} mb is more than 180000 mb of memory avaiable to thin node")
        exit(1)
    nodeN, taskN, coresN = topology(config)
    if nodeN > 1 or taskN > 1:
        print("Distributed memeory parallelization not supported")
        exit(1)
    if coresN > 80:
        print("Thin node only has 80 cores per node")
        exit(1)
    wtime = walltime(config)
    if wtime <= 30:
        partition_name = "dev_single"
    elif wtime <= 60 * 60 * 72:
        partition_name = "single"
    else:
        print(f"{sec_to_hhmmss(wtime)} is too long.")
        exit(1)
    return partition_name


config_paths = sorted_nicely(sys.argv[1:])
chain_jobs = []
for i in range(len(config_paths)):
    c = config_paths[i]
    # New chain
    if not re.search(r"chunk_[1-9]{1}\d*", c):
        chain_jobs.append([c])
    # Old chain
    else:
        chain_jobs[-1].append(c)

with open("templates/job.sh", "r") as job_file:
    job_template = job_file.read()

for chain in chain_jobs:
    jobname = re.sub(r"_chunk_\d*", rf"_chain_{len(chain)}", chain[0])
    job_file_name = Path(jobname).with_suffix(".sh").name
    jobname = Path(jobname).stem

    # Determine SLURM header on a per chain basis. Links are assumed to
    # require equal ressources
    with open(chain[0], "r") as file:
        config = json.load(file)
        mem = memory(config)
        wtime = sec_to_hhmmss(walltime(config))
        nodeN, taskN, coreN = topology(config)
        queue = partition(config)

    for config in chain:
        job_script = effify(job_template)

    for i in range(len(chain)):
        job_script += dedent(
            f"""\
            if [[ $LINK == {i} ]]; then
                ./psidm {chain[i]}
            fi
        """
        )
    slurm_dir = Path(chain[0]).parent / "slurm"
    Path.mkdir(slurm_dir, exist_ok=True)
    job_file = slurm_dir / job_file_name
    with open(job_file, "w") as file:
        file.write(job_script)
    job_file.chmod(0o744)

    # Copy submission loop script to chain job directory
    submission_file = Path("templates/submit.sh")
    submission_file_dest = slurm_dir / "submit.sh"
    copyfile(submission_file, submission_file_dest)
    submission_file_dest.chmod(0o744)
