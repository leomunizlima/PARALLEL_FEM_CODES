#include "../preconditioners.h"

int Diag_precond_CSR_setup(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I, J, neq_bef, neq, *IA, *JA;
	double *Diag, *invDiag, *AA, *f;

	AA = MatrixData->AA;
	JA = MatrixData->JA;
	IA = MatrixData->IA;
	neq = Parameters->neq;
	neq_bef = Parameters->neq_bef;
	Diag = &(MatrixData->Diag[neq_bef]);
	invDiag = &(MatrixData->invDiag[neq_bef]);
	f = &F[neq_bef];

	for (I=0; I<neq; I++){
		for (J=IA[I]; J<IA[I+1]; J++){
			if (JA[J] == I)
				Diag[I] = AA[J];	
		}
	}

	for(I=0; I<neq; I++)
		invDiag[I] = 1.0/Diag[I];	
	
	for(I=0; I<neq; I++)
		f[I] *= invDiag[I];
		
	return 0;
}



