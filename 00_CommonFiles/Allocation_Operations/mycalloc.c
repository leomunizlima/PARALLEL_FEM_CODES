#include "allocations.h"
#include <mpi.h>

void *mycalloc(char *var_name, int n, int struct_size)
{
	void *ptr;

	ptr = (void *) calloc(n,struct_size);
	
	if (ptr == NULL){
		printf("Memory allocation error in %s!\n", var_name);
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	
	return ptr;	
}



