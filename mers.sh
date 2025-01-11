#!/bin/bash

#SBATCH --job-name=Mers1
#SBATCH --mail-user=yangpe@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --account=ionides0
#SBATCH --partition=standard

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1

## 5/GB per cpu is default

#SBATCH --mem-per-cpu=1GB


## wall time is maximum time you can run

#SBATCH --time=08:00:00

### Load software modules

module load R
module list

echo "Running on $SLURM_JOB_NODELIST"
echo "Running in $(pwd)"

Rscript --vanilla mers_run.R

### sbatch script.sh

### squeue