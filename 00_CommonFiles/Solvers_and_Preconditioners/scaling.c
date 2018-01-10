#include "scaling.h"
#include "preconditioners.h"

int NO_scaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	return 0;
}


int Left_scaling_EBE(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	MPI_Request request;
	int I, J, E, nel = Parameters->nel;
	int size = NNOEL*NDOF;
	int **lm = FemStructs->lm;
	double **A = MatrixData->A;
	double *invDiag = MatrixData->invDiag;
	double invDiagEnd;	

	
	Diag_precond_EBE_setup(Parameters, MatrixData, FemStructs, 0, FemStructs->F);
	MPI_VectorUpdate(Parameters->RANK, Parameters->NPROC, Parameters, FemStructs, &invDiag[Parameters->neq_bef], 0);

	if (Parameters->RANK == Parameters->NPROC-1){
		invDiagEnd = invDiag[Parameters->NEQ-1];
		for (I=0; I<Parameters->NPROC-1;I++)
			MPI_Isend(&invDiagEnd,1,MPI_DOUBLE,I,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
	}
	else
		MPI_Recv(&invDiagEnd,1,MPI_DOUBLE,Parameters->NPROC-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	
	invDiag[Parameters->NEQ] = invDiagEnd;

	for (E=0; E<nel; E++){
		for (I=0; I<size; I++){
			for (J = 0; J<size; J++) 
				A[E][size*I+J] *= invDiag[lm[E][I]];	
		}
	}

	return 0;	
}

int Left_scaling_CSR(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	int k1, k2, I, J;
	int neq, neq_bef, neqaux_bef, neqaux_aft, *aux_bef, *aux_aft;
	double *AA, *AA_bef, *AA_aft;
	int *IA, *IA_bef, *IA_aft;
	double *invDiag;

	neq = Parameters->neq;
	neq_bef = Parameters->neq_bef;
	neqaux_bef = Parameters->neqaux_bef;
	neqaux_aft = Parameters->neqaux_aft;

	AA = MatrixData->AA;
	IA = MatrixData->IA;
	AA_bef = MatrixData->AA_bef;
	IA_bef = MatrixData->IA_bef;
	aux_bef = MatrixData->aux_bef;
	AA_aft = MatrixData->AA_aft;
	IA_aft = MatrixData->IA_aft;
	aux_aft = MatrixData->aux_aft;
	invDiag = &(MatrixData->invDiag[neq_bef]);

	Diag_precond_CSR_setup(Parameters, MatrixData, FemStructs, 0, FemStructs->F);
	
	for (I = 0; I < neqaux_bef; I++){
		k1 = IA_bef[I];
		k2 = IA_bef[I+1];
		for (J = k1; J < k2; J++)
			AA_bef[J] *= invDiag[aux_bef[I]];
	}

	for (I=0; I<neq; I++){
		for (J = IA[I]; J<IA[I+1]; J++) 
			AA[J] *= invDiag[I];	
	}
	
	for (I = 0; I < neqaux_aft; I++){
		k1 = IA_aft[I];
		k2 = IA_aft[I+1];
		for (J = k1; J < k2; J++)
			AA_aft[J] *= invDiag[aux_aft[I]];
	}


	return 0;	
}



int LeftRight_scaling_EBE(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	int I, J, E, nel = Parameters->nel;
	int neq  = Parameters->neq;
	int NEQ = Parameters->NEQ;
	int neq_bef = Parameters->neq_bef;
	int size = NNOEL*NDOF;
	int **lm = FemStructs->lm;
	double **A = MatrixData->A;
	double *Diag = MatrixData->Diag;
	double *invDiag = MatrixData->invDiag;	
	double *F = FemStructs->F;
	
	Diag_precond_EBE_setup(Parameters, MatrixData, FemStructs, 0, F);

	invDiag[NEQ] = invDiag[NEQ-1];

	for (E=0; E<nel; E++){
		for (I=0; I<size; I++){
			for (J = 0; J<size; J++) 
				A[E][size*I+J] *= sqrt(invDiag[lm[E][I]]*invDiag[lm[E][J]]);	
		}
	}
	
	for (I=neq_bef; I<neq_bef+neq; I++)
		F[I] *= sqrt(Diag[I]); //Da funcao anterior temos F divida pela diagonal Fi/dii , agora teremos Fi/sqrt(dii) pois (Fi/dii)*sqrt(di) = Fi/sqrt(dii). 

	return 0;
}


int LeftRight_scaling_CSR(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	int k1, k2, I, J;
	int neq, neq_bef, neqaux_bef, neqaux_aft, *aux_bef, *aux_aft;
	double *AA, *AA_bef, *AA_aft;
	int *IA, *IA_bef, *IA_aft;
	int *JA, *JA_bef, *JA_aft;
	double *Diag, *invDiag;
	double *F;

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
	Diag = &(MatrixData->Diag[neq_bef]);
	invDiag = &(MatrixData->invDiag[neq_bef]);
	F = &(FemStructs->F[neq_bef]);

	Diag_precond_CSR_setup(Parameters, MatrixData, FemStructs, 0, FemStructs->F);
	MPI_VectorUpdate(Parameters->RANK, Parameters->NPROC, Parameters, FemStructs, invDiag, 0);
	
	for (I = 0; I < neqaux_bef; I++){
		k1 = IA_bef[I];
		k2 = IA_bef[I+1];
		for (J = k1; J < k2; J++)
			AA_bef[J] *= sqrt(invDiag[aux_bef[I]]*invDiag[JA_bef[J]-neq_bef]);
	}
	for (I=0; I<neq; I++){
		for (J = IA[I]; J<IA[I+1]; J++) 
			AA[J] *= sqrt(invDiag[I]*invDiag[JA[J]]);	
	}
	for (I = 0; I < neqaux_aft; I++){
		k1 = IA_aft[I];
		k2 = IA_aft[I+1];
		for (J = k1; J < k2; J++)
			AA_aft[J] *= sqrt(invDiag[aux_aft[I]]*invDiag[JA_aft[J]+neq]);
	}

	for (I=0; I<neq; I++)
		F[I] *= sqrt(Diag[I]); //Da funcao anterior temos F divida pela diagonal Fi/dii , agora teremos Fi/sqrt(dii) pois (Fi/dii)*sqrt(di) = Fi/sqrt(dii). 
		

	return 0;
}

int NO_unscaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *x)
{
	return 0;
} 


int Left_unscaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *x)
{
	int I, neq = Parameters->neq; 
	double *invDiag = &MatrixData->invDiag[Parameters->neq_bef];
	
	for (I=0; I<neq; I++)
		x[I] *= sqrt(invDiag[I]);
		
	return 0;
}


