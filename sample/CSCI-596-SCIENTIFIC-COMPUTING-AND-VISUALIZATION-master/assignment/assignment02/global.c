#include "mpi.h"
#include <stdio.h>
#define NBIN 10000000
int nprocs;  /* Number of processors */
int myid;    /* My rank */

double global_sum(double partial) {
  /* Implement your own global summation here */
	doubke mydone, hisdone;
	int bitvalue, partner;
	MPI_Status status;
	mydone = partial;
	// for ( l = 0 ; l < log(double)(nprocs)/log(2.0); l++){
	for (bitvalue = 1; bitvalue < nprocs; bitvalue *= 2) {
		partner = myid ^ bitvalue;
		MPI_Send(&mydone, 1, MPI_DOUBLE, partner, bitvalue, MPI_COMM_WORLD);
		MPI_Recv(&hisdone, 1, MPI_DOUBLE, partner, bitvalue, MPI_COMM_WORLD, &status);
		mydone += hisdone;
	}
	return mydone;

}

int main(int argc, char *argv[]) {
  double partial, sum= 0.0, avg, step, x, pi, cpu1, cpu2;
  int i;
  step = 1.0 / NBIN;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  cpu1 = MPI_Wtime();
  for (i = myid; i < NBIN; i += nprocs)
  {
	  x = (i + 0.5)*step;
	  sum += 4.0 / (1.0 + x * x);
  }
  partial = sum * step;
  pi = global_sum(partial);
  cpu2 = MPI_Wtime();
  if (myid == 0) {
    printf("Value of pi = %le\n", pi);
	printf("Execution time (s) = %le\n", cpu2-cpu1);
  }

  MPI_Finalize();

  return 0;
}