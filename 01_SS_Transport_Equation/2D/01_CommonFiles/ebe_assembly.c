#include "SSTranspEquation.h"

void ebe_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double ke[3][3])
{
	double **K;

	K = MatrixData->A;

	K[E][0] = ke[0][0];	
	K[E][1] = ke[0][1];	
	K[E][2] = ke[0][2];	
	K[E][3] = ke[1][0];	
	K[E][4] = ke[1][1];	
	K[E][5] = ke[1][2];	
	K[E][6] = ke[2][0];	
	K[E][7] = ke[2][1];	
	K[E][8] = ke[2][2];	

}


