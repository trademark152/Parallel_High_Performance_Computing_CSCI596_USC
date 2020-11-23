#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:01:59
#SBATCH --output=global_pi.out
#SBATCH -A anakano_429

WORK_HOME=/project/anakano_429/tranmt/assignment3/submission_03_mt
cd $WORK_HOME

echo "##### Strong scaling #####"
mpirun -n $SLURM_NTASKS ./global_pi
mpirun -n             2 ./global_pi
mpirun -n             1 ./global_pi

echo "##### Weak scaling   #####"
mpirun -n $SLURM_NTASKS ./global_pi_iso
mpirun -n             2 ./global_pi_iso
mpirun -n             1 ./global_pi_iso
