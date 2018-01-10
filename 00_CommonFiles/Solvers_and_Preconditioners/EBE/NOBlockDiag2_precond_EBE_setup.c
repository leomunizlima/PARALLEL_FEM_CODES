#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int NOBlockDiag2_precond_EBE_setup(ParametersType *Parameters, MatrixDataType *MatrixData, NodeType *Node, ElementType *Element, int tag, double *F, int **lm)
{
	int I, J, E, I1, I2, I3;
	double **A = MatrixData->A[0];
	int neq = Parameters-> neq;
	int nel = Parameters-> nel;
	int nnodes = Parameters->nnodes;
	int **Id, *IdAux;
	double **BlockDiag, *BlockDiagAux;

	int I_diag[3][48] = {{0,1,2,3,12,13,14,15,24,25,26,27,36,37,38,39},{52,53,54,55,64,65,66,67,76,77,78,79,88,89,90,91},{104,105,106,107,116,117,118,119,128,129,130,131,140,141,142,143}};

	if (tag==1){
		IdAux = (int*) mycalloc("IdAux of 'BlockBlockDiag_precond_EBE'",4*nnodes,sizeof(int));
		Id = (int**) mycalloc("Id of 'BlockBlockDiag_precond_EBE'", nnodes,sizeof(int*));
		BlockDiagAux = (double*) mycalloc("BlockDiagAux of 'BlockBlockDiag_precond_EBE'",16*nnodes,sizeof(double));
		BlockDiag = (double**) mycalloc("BlockDiag of 'BlockBlockDiag_precond_EBE'",nnodes,sizeof(double*));

		for (I=0;I<nnodes;I++){
			BlockDiag[I] = &BlockDiagAux[16*I]; 
			Id[I] = &IdAux[4*I];
		}
		for (I=0;I<nnodes;I++){
			if (Node[I].id[0] < 0) 
				Id[I][0] = neq;
			else
				Id[I][0] = Node[I].id[0];

			if (Node[I].id[1] < 0) 
				Id[I][1] = neq; 
			else
				Id[I][1] = Node[I].id[1];
			if (Node[I].id[2] < 0) 
				Id[I][2] = neq; 
			else
				Id[I][2] = Node[I].id[2];
			if (Node[I].id[3] < 0)
				Id[I][3] = neq; 
			else
				Id[I][3] = Node[I].id[3];
		}
		MatrixData->Id = Id;
		MatrixData->BlockDiag = BlockDiag;
		MatrixData->BlockDiagAux = BlockDiagAux;
	}
	else{
		Id = MatrixData->Id;
		BlockDiag = MatrixData->BlockDiag;
		BlockDiagAux = MatrixData->BlockDiagAux;
		memset(BlockDiagAux,0,16*nnodes*sizeof(double));
	}	


	for (E=0;E<nel;E++){
		I1 = Element[E].Vertex[0];	
		I2 = Element[E].Vertex[1];	
		I3 = Element[E].Vertex[2];	
		for (I=0;I<4;I++){
			for(J=0;J<4;J++){
				BlockDiag[I1][4*I+J] += A[E][I_diag[0][4*I+J]];
				BlockDiag[I2][4*I+J] += A[E][I_diag[1][4*I+J]];
				BlockDiag[I3][4*I+J] += A[E][I_diag[2][4*I+J]];
			}
		}
	}


	return 0;
}




