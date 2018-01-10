#include "Mesh_Reordering.h"

void reading_mesh(char *Filename, int nnodes, int nel, int NDOF, double *X, double *Y, int **Type, ElementType *Element)
{
	int I, J, tag = 1; // Testing error of input
	FILE *In;		

	In = myfopen(Filename,"r");

	tag = fscanf(In, "%d", &nnodes);

	for (I=0;I<nnodes; I++){
		tag = fscanf(In,"%lf%lf",&X[I],&Y[I]); 
		for (J=0; J<NDOF; J++)
			tag = fscanf(In,"%d", &Type[I][J]);
	}


	tag = fscanf(In, "%d", &nel);

	if (NDOF==1){
		for (I=0; I<nel; I++)
			tag = fscanf(In, "%d%d%d", &(Element[I].Vertex[0]),  &(Element[I].Vertex[1]), &(Element[I].Vertex[2]));
	}
	else{
		for (I=0; I<nel; I++)
			tag = fscanf(In, "%d%d%d%d", &(Element[I].Vertex[0]),  &(Element[I].Vertex[1]), &(Element[I].Vertex[2]), &(Element[I].Type));
		
	}

	if (tag < 0){
		printf("Input data error!\n");
		exit(1);
	}

	fclose(In);

}


