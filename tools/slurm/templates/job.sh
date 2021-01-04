#!/bin/bash
#SBATCH --job-name={jobname}    
#SBATCH --nodes={nodeN}        
#SBATCH --ntasks={taskN}       
#SBATCH --cpus-per-task={coreN}
#SBATCH --mem={mem}mb           
#SBATCH --time={wtime}       
#SBATCH --output={jobname}.log 
#SBATCH --partition={queue}   

module load compiler/gnu/8.3.1
module load numlib/mkl/2020

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE
export OMP_PLACES="threads"
export OMP_PROC_BIND="close"
export OMP_DISPLAY_ENV=TRUE


cd $HOME/psidm2/build3/src


