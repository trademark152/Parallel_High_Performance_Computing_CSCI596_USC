#include "mpi.h" /* for MPI definition */
#include <stdio.h> /* for printf */

/* global variables */
int nprocs;  /* Number of processes */
int myid;    /* My rank */

/* global_sum function: Use all partial sums to return global sum */
double global_sum(double partial) {
  /* Write your hypercube algorithm here */
    /* initialization */
    double mydone, hisdone;
    int bitvalue, partner;
    MPI_Status status;

    /* pipe input to initialize output */
    mydone = partial;

    /* use this instead of log2(p) */
    for(bitvalue=1;bitvalue<nprocs;bitvalue*=2){
        partner = myid ^ bitvalue; /* XOR operation */

        // SEND my sum to partner: bitvalue is treated as communication label here
        // if receive multiple messages from multiple sources, use distinct labels as tags (bitmask or stride)
        MPI_Send(&mydone, 1, MPI_DOUBLE,partner, bitvalue, MPI_COMM_WORLD);

        // RECEIVE sum from partner: matching bitvalue is treated as communication label here   \
        // here myid and partner both acting as sender and receiver
        MPI_Recv(&hisdone, 1, MPI_DOUBLE,partner, bitvalue, MPI_COMM_WORLD, &status);
        mydone += hisdone;
    }
    return mydone;
}

/* MAIN */
int main(int argc, char *argv[]) {
  double partial, sum, avg;
  double cpu1, cpu2;
  /*mpi initialization */
  MPI_Init(&argc, &argv);/* */
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* find out my rank */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);/* number of processors  */

  partial = (double) myid;
  printf("Node %d has partial value %le\n", myid, partial);

  /* put timer to know exccution time */
  cpu1 = MPI_Wtime();
  sum = global_sum(partial);
  cpu2 = MPI_Wtime();

  /* only print out results of master run at node 0 */
  if (myid == 0) {
    avg = sum/nprocs;
    printf("Global average = %le\n", avg);
    printf("Global sum = %le\n", sum);
    printf("Execution time (s) = %le\n",cpu2-cpu1);
  }

  MPI_Finalize();
  return 0;
}
