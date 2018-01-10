#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"
/* Define preconditioner action on matrix-vector product */
/* p = A z
 * A = (LU)^{-1}
 * L * p = w
 * U * w = z
 * z out
*/
int LU_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z){
	int I, lm0, lm1, lm2;
	double *Z, z0, z1, z2; //auxiliar
	int nel = Parameters->nel;  //number of elements
	int **lm = FemStructs->lm; //array LM2, just NEQ and neq_bef<= LM < neq_bef + neq 
	double **LUe = MatrixData->LUe;
//	int neq = Parameters->neq;
	int neq_bef = Parameters->neq_bef;
	int NEQ = Parameters->NEQ;

	
	MPI_VectorUpdate(Parameters->RANK, Parameters->NPROC, Parameters, FemStructs, z, 0);
	Z = &z[-neq_bef];
	Z[NEQ] = 0; 

/*	FILE *Out;
	char Filename[200];
	sprintf(Filename,"LM_%d.dat",Parameters->RANK);
	Out=fopen(Filename, "w");
*/
	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
		lm1 = lm[I][1];
		lm2 = lm[I][2];

		z0 = Z[lm0];
		z1 = Z[lm1];
		z2 = Z[lm2];
	
/*		fprintf(Out,"LM%d->",I);
		int j;
		for (j=0;j<3;j++)
	               fprintf(Out,"%d ",lm[I][j]);
		fprintf(Out,"\n ");
		fprintf(Out,"z0=%lf z1=%lf z2=%lf\n",z0,z1,z2);
*/		

		// L * p = w
		
		z1 += - LUe[I][3] * z0; // LUe[I][4] is store as 1/LUe[I][4]
		z2 += - LUe[I][6] * z0 - LUe[I][7] * z1; // LUe[I][8] is stored as 1/LUe[8]		
		
		// U * w = z
		z2 *= LUe[I][8];
		z1 = (z1 - LUe[I][5] * z2) * LUe[I][4];
		z0 += - LUe[I][1] * z1 - LUe[I][2] * z2;

		Z[lm0] = z0;
		Z[lm1] = z1;
		Z[lm2] = z2;

		Z[NEQ] = 0;

	}

//	MPI_VectorUpdate2(Parameters->RANK, Parameters->NPROC, Parameters, FemStructs, z, 0);
/*	for (I=0;I<nel; I++){
		fprintf(Out,"LM%d->",I);
		int j;
		for (j=0;j<3;j++)
	               fprintf(Out,"%d ",lm[I][j]);
		fprintf(Out,"\n ");
	}
*/
	

//	fclose(Out);

	return 0;

}


