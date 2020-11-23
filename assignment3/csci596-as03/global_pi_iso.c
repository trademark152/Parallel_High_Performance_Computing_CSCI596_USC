#include "mpi.h"
#include <stdio.h>
#define NPERP 1000000000 // number of bin per processor
int nprocs;  /* Number of processors */
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
    double partial, sum = 0.0, avg, step, x, pi, cpu1, cpu2;
    int i;
    long long NBIN;


    /*mpi initialization */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* find out my rank */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // step size
    NBIN = (long long) NPERP*nprocs;
    step = 1.0/NBIN;

    cpu1 = MPI_Wtime();
    for(i=myid; i<NBIN; i+=nprocs){
        x = (i+0.5)*step;
        sum += 4.0/(1.0 + x*x);
    }

    // use partial as buffer value of pi for parallelization
    partial = sum*step;
    pi = global_sum(partial);
    cpu2 = MPI_Wtime();

    /* only print out results of master run at node 0 */
    if (myid == 0) {
        printf("Number of processors: %d\n", nprocs);
        printf("Value of Pi = %le\n", pi);
        printf("Execution time (s) = %le\n", cpu2 - cpu1);
    }

    MPI_Finalize();
    return 0;
}
