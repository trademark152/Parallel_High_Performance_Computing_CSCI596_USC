#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:2
#SBATCH --time=00:00:59
#SBATCH --output=pi3.out
#SBATCH --account=anakano_429

export OMP_NUM_THREADS=2

srun --cpu-bind=none -n 2 ./pi3
