#include "preconditioners.h"

int BlockDiag2_precond(ParametersType *Parameters, MatrixDataType *MatrixData, double *p, double *z)
{	
	double a,b,c,d;
	int nnodes = Parameters->nnodes;
	int neq = Parameters->neq;
	int I, **Id = MatrixData->Id;
	double **invBlockDiag = MatrixData->invBlockDiag;

	for (I=0;I<nnodes;I++){
		a = z[Id[I][0]];
		b = z[Id[I][1]];
		c = z[Id[I][2]];
		d = z[Id[I][3]];

		z[Id[I][0]] = p[Id[I][0]] + invBlockDiag[I][0]*a + invBlockDiag[I][1]*b + invBlockDiag[I][2]*c + invBlockDiag[I][3]*d;
		z[Id[I][1]] = p[Id[I][1]] + invBlockDiag[I][4]*a + invBlockDiag[I][5]*b + invBlockDiag[I][6]*c + invBlockDiag[I][7]*d;
		z[Id[I][2]] = p[Id[I][2]] + invBlockDiag[I][8]*a + invBlockDiag[I][9]*b + invBlockDiag[I][10]*c + invBlockDiag[I][11]*d;
		z[Id[I][3]] = p[Id[I][3]] + invBlockDiag[I][12]*a + invBlockDiag[I][13]*b + invBlockDiag[I][14]*c + invBlockDiag[I][15]*d;

		z[neq] = 0;
	}
	
	return 0;
}

