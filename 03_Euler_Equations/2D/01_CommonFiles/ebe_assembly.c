#include "EulerEquations.h"

#define size NNOEL*NDOF

void ebe_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*Me)[size])
{
	int i, j, J = 0;
	double **M;

	M = MatrixData->A;

	for (i = 0; i < size; i++){
		for (j = 0; j < size; j++, J++){
			M[E][J] = Me[i][j];
		}
	}

}


