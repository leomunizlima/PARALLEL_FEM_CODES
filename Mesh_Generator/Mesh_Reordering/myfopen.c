#include "Mesh_Reordering.h"

FILE *myfopen(char *FileName, char *op)
{
	FILE *InFile;

	InFile = fopen(FileName,op);
	if (InFile == NULL)
	{
		printf("File %s not found!\n", FileName);
		exit(1);
	}	
	return InFile; 
}


