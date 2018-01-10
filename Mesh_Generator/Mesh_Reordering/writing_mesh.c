#include "Mesh_Reordering.h"

void writing_mesh(int nnodes, int nel, int NDOF, int **Type, ElementType *Element, double *X, double *Y, char *Problem)
{
	int I, J;
	char Filename[200];
	FILE *Out;	

	sprintf(Filename,"Reordered_%s_%d_%d.dat",Problem,nnodes,nel);
	Out = myfopen(Filename,"w");
	fprintf(Out,"%d\n",nnodes);
	for (I=0;I<nnodes; I++){
		fprintf(Out,"%lf\t%lf\t",X[I],Y[I]); 
		for (J=0; J<NDOF; J++)
			fprintf(Out,"%d\t", Type[I][J]);
		fprintf(Out,"\n");
	}
	fprintf(Out,"%d\n",nel);

	if (NDOF==1){
		for (I=0; I<nel; I++)
			fprintf(Out, "%d\t%d\t%d\n", Element[I].Vertex[0],  Element[I].Vertex[1], Element[I].Vertex[2]);
	}
	else{
		for (I=0; I<nel; I++)
			fprintf(Out, "%d\t%d\t%d\t%d\n", Element[I].Vertex[0],  Element[I].Vertex[1], Element[I].Vertex[2], Element[I].Type);
	}

	fclose(Out);

}


