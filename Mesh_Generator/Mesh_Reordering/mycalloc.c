#include "Mesh_Reordering.h"

void *mycalloc(char *var_name, int n, int struct_size)
{
	void *ptr;

	ptr = (void *) calloc(n,struct_size);
	
	if (ptr == NULL){
		printf("Memory allocation error in %s!\n", var_name);
		exit(1);
	}
	
	return ptr;	
}



