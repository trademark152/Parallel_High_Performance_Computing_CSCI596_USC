#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:09:59
#SBATCH --output=hmd-scale.out
#SBATCH -A anakano_429

echo '8 threads'
srun -n 1 ./hmd8
echo '4 threads'
srun -n 1 ./hmd4
echo '2 threads'
srun -n 1 ./hmd2
echo '1 thread'
srun -n 1 ./hmd1