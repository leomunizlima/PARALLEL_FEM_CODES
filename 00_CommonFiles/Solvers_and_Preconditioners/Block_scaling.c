#include "scaling.h"
#include "preconditioners.h"

#define SIZE NDOF*NNOEL

int fill_Block_Ae(double [NDOF][NDOF], int, int, double *, int [SIZE][SIZE]);
int matrix_Bondary_Conditions(double [NDOF][NDOF], double [NDOF][NDOF], double [NDOF][NDOF], double [NDOF][NDOF], double [NDOF][NDOF],int, int);
int matrix_Scaling(double *, int, double *, double [NDOF][NDOF], double [NDOF][NDOF], double [NDOF][NDOF], int [SIZE][SIZE]);

int Block_scaling_EBE(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	int I, J, K, E, I1, I2, I3;
	double **A = MatrixData->A;
	int nel = Parameters-> nel;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	int I_general[SIZE][SIZE];
	double **invDiag;

	BlockDiag_precond_EBE_setup(Parameters, MatrixData, FemStructs, Parameters->iterations+1, FemStructs->F);

	invDiag = MatrixData->invBlockDiag;

	for (I=0, K=0; I<SIZE; I++)
		for (J=0; J<SIZE; J++)	
			I_general[I][J] = K++;
	
	double A11[NDOF][NDOF], A12[NDOF][NDOF], A13[NDOF][NDOF], A21[NDOF][NDOF], A22[NDOF][NDOF], A23[NDOF][NDOF], A31[NDOF][NDOF], A32[NDOF][NDOF], A33[NDOF][NDOF];
	
	for (E=0;E<nel;E++){
		I1 = Element[E].Vertex[0];	
		I2 = Element[E].Vertex[1];	
		I3 = Element[E].Vertex[2];	

		fill_Block_Ae(A11,0,0,A[E],I_general);
		fill_Block_Ae(A12,0,1,A[E],I_general);
		fill_Block_Ae(A13,0,2,A[E],I_general);
		fill_Block_Ae(A21,1,0,A[E],I_general);
		fill_Block_Ae(A22,1,1,A[E],I_general);
		fill_Block_Ae(A23,1,2,A[E],I_general);
		fill_Block_Ae(A31,2,0,A[E],I_general);
		fill_Block_Ae(A32,2,1,A[E],I_general);
		fill_Block_Ae(A33,2,2,A[E],I_general);

		/* Adjusting boundary conditions on blocks of the matrix Ae*/
		for (I=0; I<NDOF; I++){
			matrix_Bondary_Conditions(A11, A12, A13, A21, A31, Node[I1].invP_id[I], I);
			matrix_Bondary_Conditions(A22, A21, A23, A12, A32, Node[I2].invP_id[I], I);
			matrix_Bondary_Conditions(A33, A31, A32, A13, A23, Node[I3].invP_id[I], I);
		}

		/* Scaling of the Ae = Inv*Ae */	
		matrix_Scaling(A[E],0, invDiag[I1], A11, A12, A13, I_general);
		matrix_Scaling(A[E],1, invDiag[I2], A21, A22, A23, I_general);
		matrix_Scaling(A[E],2, invDiag[I3], A31, A32, A33, I_general);
	}

	return 0;
}

int fill_Block_Ae(double M[NDOF][NDOF],int I, int J, double *A, int I_g[NNOEL*NDOF][NNOEL*NDOF])
{
	int i, j;

	for (i=0; i <NDOF; i++)
		for (j=0; j<NDOF; j++)
			M[i][j] = A[I_g[NDOF*I+i][NDOF*J+j]];

	return 0;
}


int matrix_Bondary_Conditions(double M1[NDOF][NDOF], double M2[NDOF][NDOF], double M3[NDOF][NDOF], double M4[NDOF][NDOF], double M5[NDOF][NDOF], int NodeId, int tag)
{
	int i;

	if (NodeId < 0){
		
		for (i=0; i<NDOF; i++){
			M1[tag][i] = 0;
			M1[i][tag] = 0;
			M2[tag][i] = 0;
			M3[tag][i] = 0;
			M4[i][tag] = 0;
			M5[i][tag] = 0;
		}
		M1[tag][tag] = 1;
	}

	return 0;	
}


int matrix_Scaling(double *A, int tag, double *invDiag, double M1[NDOF][NDOF], double M2[NDOF][NDOF], double M3[NDOF][NDOF], int I_g[NNOEL*NDOF][NNOEL*NDOF])
{
	int i, j, k, m1, m2, m3, n, p;

	for (i=0; i<NDOF; i++){
		for (k=0; k<NDOF; k++){
			m1 = k;
			m2 = k + NDOF;
			m3 = k + 2*NDOF;
			n = NDOF*tag + i;
			A[I_g[n][m1]] = 0;
			A[I_g[n][m2]] = 0;
			A[I_g[n][m3]] = 0;
			for (j=0; j<NDOF; j++){
				p = i*NDOF+j;
				A[I_g[n][m1]] += invDiag[p]*M1[j][k];
				A[I_g[n][m2]] += invDiag[p]*M2[j][k];
				A[I_g[n][m3]] += invDiag[p]*M3[j][k];
			}
		}
	}

	return 0;

}

int NO_scaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	return 0;
}

int NO_unscaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *x)
{
	return 0;
} 


