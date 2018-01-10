#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int LU_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int nel = Parameters->nel;
	int neq_bef = Parameters->neq_bef;
	double **LUe, **A;

	if (tag == 1){
		MatrixData->LUeaux = mycalloc("LUeaux of 'LU_precond_EBE_setup'", 9*nel,sizeof(double));
		MatrixData->LUe = mycalloc("LUe of 'LU_precond_EBE_setup'",nel,sizeof(double*));
		for (I=0; I<nel; I++){
			MatrixData->LUe[I] = &(MatrixData->LUeaux[9*I]);
		}
	}
	A = MatrixData->A;
	LUe = MatrixData->LUe;

	for (I = 0; I < nel; I++){
		LUe[I][0] = 1.0;
		LUe[I][1] = A[I][1];
		LUe[I][2] = A[I][2];
		LUe[I][3] = A[I][3];
		LUe[I][4] = 1./(1.0 - A[I][1]*A[I][3]); //store 1/Lue[4] to improve LU linear solver
		LUe[I][5] = A[I][5] - A[I][2]*A[I][3];
		LUe[I][6] = A[I][6];
		LUe[I][7] = (A[I][7] - A[I][1]*A[I][6])*LUe[I][4]; 
		LUe[I][8] = 1./(1.0 - A[I][2]*A[I][6] - LUe[I][7]*LUe[I][5]); //store 1/Lue[8] to improve LU linear solver
	}

	LU_precond_EBE (Parameters, MatrixData, FemStructs, &F[neq_bef], &F[neq_bef]);

	return 0;

}

