#include "preconditioners.h"

int BlockDiag_precond(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{	
	double a,b,c,d;
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
		d = p[Id[I][3]];

		z[Id[I][0]] = invBlockDiag[I][0]*a + invBlockDiag[I][1]*b + invBlockDiag[I][2]*c + invBlockDiag[I][3]*d;
		z[Id[I][1]] = invBlockDiag[I][4]*a + invBlockDiag[I][5]*b + invBlockDiag[I][6]*c + invBlockDiag[I][7]*d;
		z[Id[I][2]] = invBlockDiag[I][8]*a + invBlockDiag[I][9]*b + invBlockDiag[I][10]*c + invBlockDiag[I][11]*d;
		z[Id[I][3]] = invBlockDiag[I][12]*a + invBlockDiag[I][13]*b + invBlockDiag[I][14]*c + invBlockDiag[I][15]*d;

		z[NEQ-neq_bef] = 0.;
	}
	
	return 0;
}

