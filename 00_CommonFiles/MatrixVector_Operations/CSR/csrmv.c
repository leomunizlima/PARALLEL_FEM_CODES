#include "../matvec.h"
#include "../../Allocation_Operations/allocations.h"

int csrmv(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType  *FemStructs, double *P, double *Q)
{
	double *AA, *AA_bef, *AA_aft, *Pa, *Qa, *P_bef, *Q_bef, *P_aft, *Q_aft;
	int i, j, k1, k2, *IA, *JA, *IA_bef, *JA_bef, *IA_aft, *JA_aft;
	int neq, neq_bef, neqaux_bef, neqaux_aft, *aux_bef, *aux_aft;
	int RANK = Parameters->RANK;
	int NPROC = Parameters->NPROC;

	neq = Parameters->neq;
	neq_bef = Parameters->neq_bef;
	neqaux_bef = Parameters->neqaux_bef;
	neqaux_aft = Parameters->neqaux_aft;

	AA = MatrixData->AA;
	IA = MatrixData->IA;
	JA = MatrixData->JA;
	AA_bef = MatrixData->AA_bef;
	IA_bef = MatrixData->IA_bef;
	JA_bef = MatrixData->JA_bef;
	aux_bef = MatrixData->aux_bef;
	AA_aft = MatrixData->AA_aft;
	IA_aft = MatrixData->IA_aft;
	JA_aft = MatrixData->JA_aft;
	aux_aft = MatrixData->aux_aft;

	MPI_VectorUpdate(RANK,NPROC, Parameters, FemStructs, &P[neq_bef], 0);

	memset(&Q[neq_bef],0,neq*sizeof(double));
	
	P_bef = &P[0];
	Q_bef = &Q[neq_bef];

	for (i = 0; i < neqaux_bef; i++){
		k1 = IA_bef[i];
		k2 = IA_bef[i+1];
		for (j = k1; j < k2; j++)
			Q_bef[aux_bef[i]] += AA_bef[j]*P_bef[JA_bef[j]];
	}

	P_aft = &P[neq + neq_bef];
	Q_aft = &Q[neq_bef];

	for (i = 0; i < neqaux_aft; i++){
		k1 = IA_aft[i];
		k2 = IA_aft[i+1];
		for (j = k1; j < k2; j++)
			Q_aft[aux_aft[i]] += AA_aft[j]*P_aft[JA_aft[j]];
	}

	
	Pa = &P[neq_bef];
	Qa = &Q[neq_bef];
	
	for (i = 0; i < neq; i++){
		k1 = IA[i];
		k2 = IA[i+1];
		for (j = k1; j < k2; j++)
			Qa[i] += AA[j]*Pa[JA[j]];
	}

	return 0;
}








