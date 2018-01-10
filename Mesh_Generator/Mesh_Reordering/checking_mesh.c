#include "Mesh_Reordering.h"

void checking_mesh(int *nnodes_out, int nel, double *X, double *Y, int **Type, ElementType *Element)
{
	int I, J, I1, I2, I3;
	int nnodes = *nnodes_out;
	int *check = mycalloc("check of 'main'",nnodes,sizeof(int));

	for (I=0;I<nel;I++){
		I1 = Element[I].Vertex[0];
		I2 = Element[I].Vertex[1];
		I3 = Element[I].Vertex[2];

		check[I1] = 1;
		check[I2] = 1;
		check[I3] = 1;

	}

	I=0; 
	while(I<nnodes){
		if (!check[I]){
			J = I;
			while(J<nnodes-1){
				X[J] = X[J+1];
				Y[J] = Y[J+1];
				Type[J] = Type[J+1];
				J++;
			}
			nnodes--;
			for (J=0;J<nel;J++){
				if (Element[J].Vertex[0]>=I)
					Element[J].Vertex[0]--;
				if (Element[J].Vertex[1]>=I)
					Element[J].Vertex[1]--;
				if (Element[J].Vertex[2]>=I)
					Element[J].Vertex[2]--;
			}		
				
		}
		I++;
	}

	*nnodes_out = nnodes;
	free(check);

}


