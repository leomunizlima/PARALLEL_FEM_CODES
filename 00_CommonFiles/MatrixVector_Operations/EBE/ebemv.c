#include "../matvec.h"
#include "../../Allocation_Operations/allocations.h"

int ebemv(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType  *FemStructs, double *P, double *Q)
{

	int I, i1, i2, i3, NEQ, nel, neq_bef, RANK, NPROC;
	double p1, p2, p3;
	double **A, *Pa, *Qa; // Pa and Qa are used to reduce cache miss
	int **lm, *Perm, *invPerm;

	NEQ = Parameters->NEQ;
	nel = Parameters->nel;
	neq_bef = Parameters->neq_bef;
	lm = FemStructs->lm;
	A = MatrixData->A;
	Perm = MatrixData->Perm;
	invPerm = MatrixData->invPerm;
	RANK = Parameters->RANK;
	NPROC = Parameters->NPROC;

/*	char filename[200];
	FILE *Out;

	sprintf(filename,"RANK_%d.dat",RANK);
	Out = fopen(filename,"w");

	fprintf(Out,"P_a:\n");
	for (I=0; I<NEQ; I++)
		fprintf(Out,"Pa[%d]=%lf;\n", I,P[I]);
*/
	MPI_VectorUpdate(RANK,NPROC, Parameters, FemStructs, &P[neq_bef], 0);
	
	Pa = (double *) mycalloc("Pa of 'mvebe'", NEQ+1, sizeof(double));
	Qa = (double *) mycalloc("Qa of 'mvebe'", NEQ+1, sizeof(double));
 
	for (I=0; I < NEQ; I++)
		Pa[I] = P[Perm[I]];

	for(I = 0; I < nel; I++){

		i1 = invPerm[lm[I][0]];
		i2 = invPerm[lm[I][1]];
		i3 = invPerm[lm[I][2]];

		p1 = Pa[i1];		
		p2 = Pa[i2];		
		p3 = Pa[i3];		
			
		Qa[i1] += A[I][0]*p1 + A[I][1]*p2 + A[I][2]*p3;
		Qa[i2] += A[I][3]*p1 + A[I][4]*p2 + A[I][5]*p3;
		Qa[i3] += A[I][6]*p1 + A[I][7]*p2 + A[I][8]*p3;
				
		Qa[NEQ] = 0;
	}

	for (I = 0; I < NEQ; I++)
		Q[Perm[I]] = Qa[I];

	free(Pa);
	free(Qa);	

	return 0;
}


