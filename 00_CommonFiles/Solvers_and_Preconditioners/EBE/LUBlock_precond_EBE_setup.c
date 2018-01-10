#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int LUBlock_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I, J, j, K;
	int nel = Parameters->nel;
	int size = NDOF*NNOEL;
	int size2 = size*size;
	double A[NDOF*NNOEL][NDOF*NNOEL];

	if (tag == 1){
		MatrixData->LUeaux = mycalloc("LUeaux of 'LUBlock_precond_EBE_setup'", size2*nel,sizeof(double));
		MatrixData->LUe = mycalloc("LUe of 'LUBlock_precond_EBE_setup'",nel,sizeof(double*));
		for (I=0; I<nel; I++){
			MatrixData->LUe[I] = &(MatrixData->LUeaux[size2*I]);
		}
	}

	double m;
	for (K=0; K<nel; K++)
	{
		memcpy (&A[0][0], MatrixData->A[K], size2*sizeof(double));
		for (I=0; I<NDOF; I++){
			for (J=0; J<NDOF; J++){
				A[I][J] = 0;
				A[I+NDOF][J+NDOF] = 0;
				A[I+2*NDOF][J+2*NDOF] = 0;
			}
				
		}
		for (I=0; I<NDOF; I++){
			A[I][I] = 1.0;
			A[I+NDOF][I+NDOF]  = 1.0;	
			A[I+2*NDOF][I+2*NDOF]  = 1.0;	
		}
			
		for (I=0; I<NDOF; I++){
			for (J=NDOF; J<size; J++){
				m = A[J][I];
				for (j=I+1; j<size; j++){		
					A[J][j] += -m*A[I][j];
				}
			}
		}
		
		for (I=NDOF; I<size-1; I++){
			for (J=NDOF; J<size; J++){
				m = A[J][I]/A[I][I];
				for (j=I+1; j<size; j++){		
					A[J][j] += -m*A[I][j];
				}
				A[J][I] = m;
			}
		}

		for (I=0; I<size; I++)
			A[I][I] = 1.0/A[I][I];

		
		memcpy (MatrixData->LUe[K], &A[0][0], size2*sizeof(double));
	}	


/*	for (I = 0; I < nel; I++){
//		LUe[I][0] = 1.0/sqrt(1.0-A[I][48]*A[I][4]);
//		LUe[I][1] = 1.0/sqrt(1.0-A[I][61]*A[I][17]);
//		LUe[I][2] = 1.0/sqrt(1.0-A[I][74]*A[I][30]);
//		LUe[I][3] = 1.0/sqrt(1.0-A[I][87]*A[I][43]);
		LUe[I][0] = 1.0/(1.0 - A[I][48]*A[I][4]);
		LUe[I][1] = 1.0/(1.0 - A[I][61]*A[I][17]);
		LUe[I][2] = 1.0/(1.0 - A[I][74]*A[I][30]);
		LUe[I][3] = 1.0/(1.0 - A[I][87]*A[I][43]);

		LUe[I][4] = A[I][56] - A[I][48]*A[I][8];
		LUe[I][5] = A[I][69] - A[I][61]*A[I][21];
		LUe[I][6] = A[I][82] - A[I][74]*A[I][34];
		LUe[I][7] = A[I][95] - A[I][87]*A[I][47];

		LUe[I][8] = 1.0 - A[I][96]*A[I][4]*LUe[I][0];
		LUe[I][9] = 1.0 - A[I][109]*A[I][17]*LUe[I][1];
		LUe[I][10] = 1.0 - A[I][122]*A[I][30]*LUe[I][2];
		LUe[I][11] = 1.0 - A[I][135]*A[I][43]*LUe[I][3];

		LUe[I][12] = 1./(1.0 - A[I][96]*A[I][8] - LUe[I][8]); //store 1/Lue[8] to improve LU linear solver
		LUe[I][13] = 1./(1.0 - A[I][109]*A[I][21] - LUe[I][9]); //store 1/Lue[8] to improve LU linear solver
		LUe[I][14] = 1./(1.0 - A[I][122]*A[I][34] - LUe[I][10]); //store 1/Lue[8] to improve LU linear solver
		LUe[I][15] = 1./(1.0 - A[I][135]*A[I][47] - LUe[I][11]); //store 1/Lue[8] to improve LU linear solver
	}
*/
		

	LUBlock_precond_EBE (Parameters, MatrixData, FemStructs, F, F);

	return 0;

}

