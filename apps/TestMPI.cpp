/*The Parallel Hello World Program*/
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

int main(int argc, char **argv)
{
  int node;
  int provided;
  provided = MPI::Init_thread(MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD, &node);
  printf("Hello World from Node %d\n",node);
  printf("Supports level %d of %d %d %d %d\n", provided, MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE);
  MPI::Finalize();
  return 0;
}
 
