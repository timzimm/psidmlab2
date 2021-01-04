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

cd $HOME/psidm2/build2/src


