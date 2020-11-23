#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:04:59
#SBATCH --output=pmd_split.out
#SBATCH -A anakano_429

WORK_HOME=/scratch2/tranmt/csci_589/assignment4
cd $WORK_HOME

mpirun -n $SLURM_NTASKS ./pmd_split
