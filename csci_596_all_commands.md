# Week 1
## Run a program
  gcc -o lmd_sqrt lmd_sqrt.c -lm
  lmd_sqrt < lmd.in


# Week 2
## Discovery cluster
  pwd: directory
  ls: file listings
  ls -a: all files (hidden)
  vim .bashrc: open text editor to edit bashrc, to exit ":q", to save :w"
  less .bashrc: view file
  echo $0: view shell type
  which mpicc
  which mpirun
  myaccount: all accounts to charge
  source .bashrc: specify source
  mpicc -o mpi_simple mpi_simple.c: compile mpi program
  sbatch mpi_simple.sl
  myquota

## bashrc edit: to set up standard software environment
    # CUDA
    module purge
    module load usc
    module load cuda/10.1.243

## run jobs on slurm batch system: typical slurm files
    #!/bin/bash  -> type of program to run
    #SBATCH --ntasks-per-node=2  -> request 2 processors per node
    #SBATCH --nodes=1  -> require 1 computing node
    #SBATCH --time=00:00:10
    #SBATCH --output=mpi_simple.out -> output file name
    #SBATCH -A anakano_429 -> account name
    WORK_HOME=/scratch2/tranmt/csci_596
    cd $WORK_HOME
    mpirun -n $SLURM_NTASKS ./mpi_simple

## interactive job
    salloc -n 2 -t 20 -> reserve 2 processors for 20 minutes
    mpirun -n 2 ./mpi_simple
    less /proc/cpuinfo -> find what kind of node you get
    exit -> return to fron-end discovery


## symbolic link to project folder
    ln -s /scratch2/tranmt/csci_596 cs596
    ls -lt
    cd cs596
    pwd -P

## transfer files
    sftp tranmt@discovery.usc.edu
    cd "project_file_directoryxit"
    put "file"
    get "file"

## MPI program 1
    #include "mpi.h"  -> MPI definition (available in cluster)
    #include <stdio.h>  -> for printf
    int main(int argc, char *argv[]) {
    MPI_Status status;
    int myid;
    int n;
    MPI_Init(&argc, &argv); -> register your running program
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); -> return MPI rank (process ID). MPI_COMM_WORLD contains all the running processes
    if (myid == 0) {
      n = 777;
      MPI_Send(&n, 1, MPI_INT, 1, 10, MPI_COMM_WORLD); -> &n (starting address of data), 1(# of data items), MPI_INT (data_type), these three forming a data triplet; 1 is source/destination rank (to/from whom), 10 (tag/label of message subject), MPI_COMM_WORLD (communicator)
    }
    else {
      MPI_Recv(&n, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status); -> &status is to
      printf("n = %d\n", n);
    }
    MPI_Finalize(); -> de-register your running program
    return 0;
    }

## MPI global operation
    int l_v, g_s; // local variable & global sum
    l_v = myid; // myid is my MPI rank
    MPI_Allreduce(&l_v, &g_s, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

## MPI barrier
    MPI_Barrier(MPI_Comm communicator)

## MPI slurm files mpi_simple.sl
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=2
    #SBATCH --time=00:00:10
    #SBATCH --output=mpi_simple.out
    mpirun -n $SLURM_NTASKS ./mpi_simple

## MPI program
    mpicc -o mpi_simple mpi_simple.c -> compile an MPI program
    mpirun -n 2 mpi_simple -> execute an MPI program
    which mpicc -> to find absolute path to mpicc
    more /proc/cpuinfo -> to find processor information

## MPI on HPC
    sbatch mpi_simple.sl
    squeue -user tranmt
    scancel jobID
    more mpi_simple.out -> check output
    salloc -n 2 -t 20 -> run interactively to debug: 2 processors for 20 minutes
    mpirun -n 2 ./mpi_simple -> 2 means 2 instances of the same program mpi_simple will be spawned simultaneously
    exit
    ln -s /scratch2/tranmt/csci_596 cs596 -> symbolic link
    sftp tranmt@discovery.usc.edu -> log in to transfer file (alternative to fileZilla)
    cd cs596
    put md.* -> transfer from local to remote computer
    ls -> check if the file has been transfered
    get md.* ->> transfer file from remote to local

# Week 3
## mpi_comm.c:
    ping hpc-transfer.usc.edu -> calculate latency
    traceroute www.u-tokyo-ac.jp -> trace route of data routing

## Asynchronous message passing
    MPI_Isend()
    MPI_Irecv()
    MPI_Wait
    MPI_Waitall
    MPI_Waitany
    MPI_Test

## parallel MD
    cat remd22.f | grep MPI | less -> | is a pipe, grep is a pattern matching command

# Week 4
## Local to GIT
git init -> initialize
git add README.md or git add . -> add all files to staging
git commit -m "first commit" -> commit with comment
git branch -M master -> add mastering branch
git remote add origin https://github.com/trademark152/repository_name.git -> add address of GIT remote, alias as origin
git push -u origin master -> push data to GIT repo

## Remote to GIT
1) initialization

git --help
~ $ git config –global user.name “XXX”
~ $ git config –global user.email “your_ID@usc.edu”
~ $ git config –global core.editor “vim”

2) Create local folder on discovery node
cd to repo_folder
~/repo_folder $ git init -> Create an empty Git repository
~/repo_folder $ git add * -> Stage all files in the directory to be tracked by Git
~/repo_folder $ git commit -> Record changes to the repository
~/repo_folder $ git remote origin https://github.com/trademark152?repo_name.git
~/repo_folder $ git push origin master
~/ $ git clone https://github.com/trademark152/repo_name.git -> clone repo
~/ $ git pull origin master -> retrieve updated commits

## GIT commands
