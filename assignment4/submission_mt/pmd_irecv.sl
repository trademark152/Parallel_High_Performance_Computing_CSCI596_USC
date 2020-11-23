#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:01:59
#SBATCH --output=pmd_irecv.out
#SBATCH -A anakano_429

counter=0
while [ $counter -lt 3 ]; do
  echo "***** Asynchronous *****"
  mpirun -n $SLURM_NTASKS ./pmd_irecv
  echo "***** Synchronous *****"
  mpirun -n $SLURM_NTASKS ./pmd
  let counter+=1
done
