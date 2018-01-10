#include "preconditioners.h"

int NOBlockDiag2_precond(ParametersType *Parameters, MatrixDataType *MatrixData, double *p, double *z)
{	
	double a,b,c,d;
	int nnodes = Parameters->nnodes;
	int neq = Parameters->neq;
	int I, **Id = MatrixData->Id;
	double **BlockDiag = MatrixData->BlockDiag;

	for (I=0;I<nnodes;I++){
		a = p[Id[I][0]];
		b = p[Id[I][1]];
		c = p[Id[I][2]];
		d = p[Id[I][3]];

		z[Id[I][0]] += BlockDiag[I][0]*a + BlockDiag[I][1]*b + BlockDiag[I][2]*c + BlockDiag[I][3]*d;
		z[Id[I][1]] += BlockDiag[I][4]*a + BlockDiag[I][5]*b + BlockDiag[I][6]*c + BlockDiag[I][7]*d;
		z[Id[I][2]] += BlockDiag[I][8]*a + BlockDiag[I][9]*b + BlockDiag[I][10]*c + BlockDiag[I][11]*d;
		z[Id[I][3]] += BlockDiag[I][12]*a + BlockDiag[I][13]*b + BlockDiag[I][14]*c + BlockDiag[I][15]*d;

		z[neq] = 0;
	}
	
	return 0;
}


