#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=100
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=standard
##BATCH --reservation=shortqos
#SBATCH --chdir=../run

# Setup the job environment (this module needs to be loaded before any other modules)
module load epcc-job-env

# Set the number of threads to 1
export OMP_NUM_THREADS=1

# Set scratch directory
export TMPDIR=/work/n01/n01/ryapat01/SCRATCH

# run the job 
srun --distribution=block:block --hint=nomultithread ./mitgcmuv

