#include "preconditioners.h"

int BlockDiagDOF3_precond(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{	
	double a,b,c;
	int nnodes = Parameters->nnodes;
	int NEQ = Parameters->NEQ;
	int neq_bef = Parameters->neq_bef;
	int I, **Id = MatrixData->Id;
	double **invBlockDiag = MatrixData->invBlockDiag;

	p[NEQ-neq_bef] = 0;

	for (I=0;I<nnodes;I++){
		a = p[Id[I][0]];
		b = p[Id[I][1]];
		c = p[Id[I][2]];

		z[Id[I][0]] = invBlockDiag[I][0]*a + invBlockDiag[I][1]*b + invBlockDiag[I][2]*c;
		z[Id[I][1]] = invBlockDiag[I][3]*a + invBlockDiag[I][4]*b + invBlockDiag[I][5]*c;
		z[Id[I][2]] = invBlockDiag[I][6]*a + invBlockDiag[I][7]*b + invBlockDiag[I][8]*c;

		z[NEQ-neq_bef] = 0.;
	}
	
	return 0;
}

