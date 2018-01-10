#include "TranspEquation.h"

int setzeros(ParametersType *Parameters, MatrixDataType *MatrixData)
{
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
		int size = NDOF*NNOEL;
		int size2 = size*size;

		memset(MatrixData->Aaux,0,Parameters->nel*size2*sizeof(double));
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		memset(MatrixData->AA,0,(Parameters->nnzero + 1)*sizeof(double));
		memset(MatrixData->AA_bef,0,(Parameters->nnzero_bef + 1)*sizeof(double));
		memset(MatrixData->AA_aft,0,(Parameters->nnzero_aft + 1)*sizeof(double));
	}
	return 0;
}

