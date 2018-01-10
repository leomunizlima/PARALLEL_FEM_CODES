#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int DOF3_matrix_Inverse(double *, double *);

int BlockDiagDOF3_precond_EBE_setup(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I, J, E, I1, I2, I3;
	double **A = MatrixData->A;
	int NEQ = Parameters-> NEQ;
	int neq_bef = Parameters->neq_bef;
	int nel = Parameters-> nel;
	int nnodes = Parameters->nnodes;
	int **Id, *IdAux;
	double **BlockDiag, *BlockDiagAux, **invBlockDiag, *invBlockDiagAux;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	int I_diag[3][9] = {{0,1,2,9,10,11,18,19,20},{30,31,32,39,40,41,48,49,50},{60,61,62,69,70,71,78,79,80}};

	if (tag==1){	
		IdAux = (int*) mycalloc("IdAux of 'BlockBlockDiag_precond_EBE'",3*nnodes,sizeof(int));
		Id = (int**) mycalloc("Id of 'BlockBlockDiag_precond_EBE'", nnodes,sizeof(int*));
		BlockDiagAux = (double*) mycalloc("BlockDiagAux of 'BlockBlockDiag_precond_EBE'",9*nnodes,sizeof(double));
		BlockDiag = (double**) mycalloc("BlockDiag of 'BlockBlockDiag_precond_EBE'",nnodes,sizeof(double*));
		invBlockDiagAux = (double*) mycalloc("invBlockDiagAux of 'BlockBlockDiag_precond_EBE'",9*nnodes,sizeof(double));
		invBlockDiag = (double**) mycalloc("invBlockDiag of 'BlockBlockDiag_precond_EBE'",nnodes,sizeof(double*));
		for (I=0;I<nnodes;I++){
			BlockDiag[I] = &BlockDiagAux[9*I]; 
			invBlockDiag[I] = &invBlockDiagAux[9*I]; 
			Id[I] = &IdAux[3*I];
		}
		for (I=0;I<nnodes;I++){
			if (Node[I].Type[0] == 0) 
				Id[I][0] = NEQ-neq_bef;
			else
				Id[I][0] = Node[I].invP_id[0];

			if (Node[I].Type[1] == 0) 
				Id[I][1] = NEQ-neq_bef; 
			else
				Id[I][1] = Node[I].invP_id[1];
			if (Node[I].Type[2] == 0) 
				Id[I][2] = NEQ-neq_bef; 
			else
				Id[I][2] = Node[I].invP_id[2];
		}
		MatrixData->BlockDiag = BlockDiag;
		MatrixData->BlockDiagAux = BlockDiagAux;
		MatrixData->invBlockDiag = invBlockDiag;
		MatrixData->Id = Id;
	}
	else{
		Id = MatrixData->Id;
		BlockDiag = MatrixData->BlockDiag;
		BlockDiagAux = MatrixData->BlockDiagAux;
		memset(BlockDiagAux,0,9*nnodes*sizeof(double));
		invBlockDiag = MatrixData->invBlockDiag;
	}	

	for (E=0;E<nel;E++){
		I1 = Element[E].Vertex[0];	
		I2 = Element[E].Vertex[1];	
		I3 = Element[E].Vertex[2];	
		for (I=0;I<3;I++){
			for(J=0;J<3;J++){
				BlockDiag[I1][3*I+J] += A[E][I_diag[0][3*I+J]];
				BlockDiag[I2][3*I+J] += A[E][I_diag[1][3*I+J]];
				BlockDiag[I3][3*I+J] += A[E][I_diag[2][3*I+J]];
			}
		}
	}

	for (I=0;I<nnodes;I++)
		DOF3_matrix_Inverse(BlockDiag[I],invBlockDiag[I]);

	/* Preconditioning of F */
	BlockDiagDOF3_precond(Parameters, MatrixData, FemStructs,&F[neq_bef], &F[neq_bef]);
				

	return 0;
}

int DOF3_matrix_Inverse(double *BlockDiag, double *invBlockDiag)
{
	double a, b,c,d,e,f,g,h,i,det;

	/*This function calculates the inverse of a 3 x 3 matrix 
	
	   |a b c|	
	A =|d e f| 	
	   |g h i|	*/
	  

	a = BlockDiag[0];
	b = BlockDiag[1];
	c = BlockDiag[2];
	d = BlockDiag[3];
	e = BlockDiag[4];
	f = BlockDiag[5];
	g = BlockDiag[6];
	h = BlockDiag[7];
	i = BlockDiag[8];
	
	/* Determinant calculations */
	det = a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g;

//	if (fabs(det)<1e-16) printf("OPS!!!! Determante NULO!!!\n");

	double invdet = 1.0/det;

	/* Inverse of matrix calculations */
	invBlockDiag[0] = invdet*(e*i-f*h);		 	
	invBlockDiag[1] = invdet*(c*h-b*i);		 	
	invBlockDiag[2] = invdet*(b*f-c*e);		 	
	invBlockDiag[3] = invdet*(f*g-d*i);		 	
	invBlockDiag[4] = invdet*(a*i-c*g);		 	
	invBlockDiag[5] = invdet*(c*d-a*f);		 	
	invBlockDiag[6] = invdet*(d*h-e*g);		 	
	invBlockDiag[7] = invdet*(b*g-a*h);		 	
	invBlockDiag[8] = invdet*(a*e-b*d);		 	

	return 0;
}		



