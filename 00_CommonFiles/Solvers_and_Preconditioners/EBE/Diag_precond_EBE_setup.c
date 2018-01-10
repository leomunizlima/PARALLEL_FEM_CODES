#include "../preconditioners.h"

int Diag_precond_EBE_setup(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int E, I, K, nel, neq, NEQ, neq_bef, size=NDOF*NNOEL; 
	int **lm = FemStructs->lm;
	double *Diag, *invDiag, **A;

	neq = Parameters-> neq;
	NEQ = Parameters-> NEQ;
	neq_bef = Parameters->neq_bef;
	nel = Parameters-> nel;
	A = MatrixData->A;
	Diag = MatrixData->Diag;
	invDiag = MatrixData->invDiag;

	for (E=0; E<nel; E++){
		for (I=0, K=0; I<size; I++){
			Diag[lm[E][I]] += A[E][K];
			K += size+1;
		}
		Diag[NEQ] = 0;
	}

	for(I=neq_bef; I<neq+neq_bef; I++)
		invDiag[I] = 1.0/Diag[I];	

	for(I=neq_bef; I<neq+neq_bef; I++)
		F[I] *= invDiag[I];
			
	return 0;
}


