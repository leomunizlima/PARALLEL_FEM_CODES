#include "preconditioners.h"

int Diag_precond(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{		
	int I, neq = Parameters->neq;
	int neq_bef = Parameters->neq_bef;
	double *invDiag = &(MatrixData->invDiag[neq_bef]);


	for (I=0; I<neq; I++)
		z[I] = p[I]*invDiag[I];


	return 0;
}
	
