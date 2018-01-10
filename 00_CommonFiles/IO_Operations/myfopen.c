#include "io.h"
#include <mpi.h>

FILE *myfopen(char *FileName, char *op)
{
	FILE *InFile;

	InFile = fopen(FileName,op);
	if (InFile == NULL)
	{
		printf("File %s not found!\n", FileName);
		MPI_Abort(MPI_COMM_WORLD,1);
	}	
	return InFile; 
}


