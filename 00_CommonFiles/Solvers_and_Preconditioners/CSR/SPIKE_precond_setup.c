#include "../preconditioners.h"
#include "../pardiso.h"
#include "../pardiso.c"
#include "../../Allocation_Operations/allocations.h"
#include "mpi.h"
#include <limits.h>

void PARALLEL_REDUCED_system (double* V, double* W, int k, MAT* S, int RANK, int tag);

int SPIKE_precond_setup(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int k, I, J; 
	int n = Parameters->neq;
	int n_bef = Parameters->neq_bef;
	int n_aft = Parameters->NEQ - n_bef - n;
	int naux_aft = Parameters->neqaux_aft;
	int nz = Parameters->nnzero;
	int RANK = Parameters->RANK;
	int NPROC = Parameters->NPROC;
	double **fB=NULL;
	double **fC=NULL;
	double *fBaux=NULL;
	double *fCaux=NULL;
	double *V, *W, *WT;
	MAT *A, *B, *C, *S;
	PardisoVariablesType *PardisoVariables;
	void** A_pt; 
	int* A_iparm; 
	double* A_dparm; 
	void** S_pt; 
	int* S_iparm; 
	double* S_dparm;
	MPI_Request request;
		
	if (tag == 1){	 
		int N_bef;
		int N_aft;
		int temp;

		if (n_bef==0)
			temp = INT_MAX;
		else
			temp = n_bef;

		MPI_Allreduce(&temp,&N_bef,1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

		if (n_aft==0)
			temp = INT_MAX;
		else
			temp = n_aft;

		MPI_Allreduce(&temp,&N_aft,1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

		k = atoi(&(Parameters->Preconditioner[5])); // In SPIKE, k represents the size of coupling blocks
		k = min(k,N_bef);
		k = min(k,N_aft);
		if (k != atoi(&(Parameters->Preconditioner[5])))
			sprintf(&(Parameters->Preconditioner[5]),"%d",k);

		printf ("=========> k=%d (RANK=%d) <=====\n",k, RANK);

		A = mycalloc("A of 'SPIKE_precond_setup'",1,sizeof(MAT));
		S = mycalloc("S of 'SPIKE_precond_setup'",1,sizeof(MAT));
		C = mycalloc("C of 'SPIKE_precond_setup'",1,sizeof(MAT));
		B = mycalloc("B of 'SPIKE_precond_setup'",1,sizeof(MAT));
		
		A->AA  = MatrixData->AA;
		A->D   = &(MatrixData->Diag[n_bef]);
		A->JA  = MatrixData->JA;
		A->IA  = MatrixData->IA;
		A->n   = n;
		A->m   = n;
		A->nz  = nz;

		if (RANK != 0){
			fCaux = mycalloc("fCaux of 'SPIKE_precond_setup'", n*k, sizeof(double)); 
			fC    = mycalloc("fC of 'SPIKE_precond_setup'", n, sizeof(double*)); 
		
			for (I = 0; I<n; I++)
				fC[I] = &fCaux[I*k];

		}
		MatrixData->fullCaux = fCaux;
		MatrixData->fullC = fC;

		if (RANK != NPROC-1){	
			fBaux = mycalloc("fBaux of 'SPIKE_precond_setup'", n*k, sizeof(double)); 
			fB    = mycalloc("fB of 'SPIKE_precond_setup'", n, sizeof(double*)); 
		
			for (I = 0; I<n; I++)
				fB[I] = &fBaux[I*k];

	
		}
		MatrixData->fullBaux = fBaux;
		MatrixData->fullB = fB;
		
		MatrixData->A_SPIKE = A;
		MatrixData->S_SPIKE = S;
		MatrixData->C_SPIKE = C;
		MatrixData->B_SPIKE = B;
		Parameters->k_SPIKE_block_size = k;
		PardisoVariables =  mycalloc("PardisoVariables of 'SPIKE_precond_setup'",1,sizeof(PardisoVariablesType));
		Parameters->PardisoVariables = PardisoVariables;
	}

//	if (tag-1%Parameters->NonLinearMaxIter!=0) return 0;	
	
	k = Parameters->k_SPIKE_block_size;
	fC = MatrixData->fullC;
	fB = MatrixData->fullB;
	A = MatrixData->A_SPIKE;
	S = MatrixData->S_SPIKE;
	PardisoVariables = Parameters->PardisoVariables;
	A_pt 	= PardisoVariables->A_pt; 
	A_iparm	= PardisoVariables->A_iparm; 
	A_dparm	= PardisoVariables->A_dparm; 
	S_pt 	= PardisoVariables->S_pt; 
	S_iparm	= PardisoVariables->S_iparm; 
	S_dparm	= PardisoVariables->S_dparm; 

	PARDISO_release_memory (A, A_pt, A_iparm, A_dparm, 1);
	PARDISO_release_memory (S, S_pt, S_iparm, S_dparm, 1);
				

	if (RANK != 0){
		double 	*AA_bef = MatrixData->AA_bef;
		int	*JA_bef = MatrixData->JA_bef;
		int	*IA_bef = MatrixData->IA_bef;
		int 	offset = max(n_bef - k,0);
		int 	nz_C = 0;

		for (I=0; I<k; I++)
			for (J = IA_bef[I]; J<IA_bef[I+1]; J++)			
				if(JA_bef[J] >= offset){ 
					fC[I][JA_bef[J]-offset] = AA_bef[J];
					nz_C++;
				}
		
		if (tag==1){
			C->AA 	= mycalloc("C->AA of 'SPIKE_precond_setup'",nz_C,sizeof(double));
			C->JA 	= mycalloc("C->AA of 'SPIKE_precond_setup'",nz_C,sizeof(int));
			C->IA 	= mycalloc("C->AA of 'SPIKE_precond_setup'",k+1,sizeof(int));
			C->D	= mycalloc("C->D of 'SPIKE_precond_setup'",k,sizeof(double));
			C->n 	= k;
			C->m	= k;
			C->nz	= nz_C;
		}
		else
			C = MatrixData->C_SPIKE;

		nz_C = 0;
		for (I=0; I<k; I++)
			for (J = IA_bef[I]; J<IA_bef[I+1]; J++)			
				if(JA_bef[J] >= offset){ 
					C->AA[nz_C] = AA_bef[J];
					C->JA[nz_C] = JA_bef[J]-offset;	
					C->IA[I] = nz_C;
					if (I==C->JA[nz_C])
						C->D[I] = AA_bef[J];	
					nz_C++;
				}
		C->IA[k] = nz_C;

	}	

	if (RANK != NPROC-1){
		double	*AA_aft = MatrixData->AA_aft;
		int	*JA_aft = MatrixData->JA_aft;
		int 	*IA_aft = MatrixData->IA_aft;
		int 	offset = max(naux_aft-k,0);
		int	nz_B = 0;

		for (I=offset; I<naux_aft; I++)
			for (J = IA_aft[I]; J<IA_aft[I+1]; J++)			
				if(JA_aft[J] < k){  
					fB[I-offset][JA_aft[J]] = AA_aft[J];
					nz_B++;
				}
		if (tag==1){
			B->AA 	= mycalloc("B->AA of 'SPIKE_precond_setup'",nz_B,sizeof(double));
			B->JA 	= mycalloc("B->AA of 'SPIKE_precond_setup'",nz_B,sizeof(int));
			B->IA 	= mycalloc("B->AA of 'SPIKE_precond_setup'",k+1,sizeof(int));
			B->D 	= mycalloc("B->AA of 'SPIKE_precond_setup'",k,sizeof(double));
			B->n 	= k;
			B->m	= k;
			B->nz	= nz_B;
		}
		else
			B = MatrixData->B_SPIKE;

		nz_B = 0;
		for (I=offset ; I<naux_aft; I++)
			for (J = IA_aft[I]; J<IA_aft[I+1]; J++)			
				if(JA_aft[J] < k){  
					B->AA[nz_B] = AA_aft[J];
					B->JA[nz_B] = JA_aft[J];	
					B->IA[I-offset] = nz_B;
					if (I-offset==B->JA[nz_B])
						B->D[I-offset] = AA_aft[J];
					nz_B++;
				}
		B->IA[k] = nz_B;

	}

	
	if (RANK == 0)
	{
		V  = mycalloc ("V of SPIKE_precond_setup'",k*k ,sizeof(double));
		WT = mycalloc ("WT of SPIKE_precond_setup'",k*k ,sizeof(double));

		PARDISO_bottom_tips (A,&fB[0][0],V,k);

		MPI_Recv(WT,k*k,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		PARALLEL_REDUCED_system (V,WT,k,S,RANK, tag);

		free(V);
		free(WT);

		PARDISO_numerical_factorization (A,A_pt,A_iparm,A_dparm,1);
		PARDISO_numerical_factorization (S,S_pt,S_iparm,S_dparm,1);
	}
	else if (RANK == NPROC - 1)
	{
		W  = mycalloc ("W of SPIKE_precond_setup'",k*k ,sizeof(double));

		PARDISO_top_tips  (A,&fC[0][0],W,k);

		MPI_Isend(W,k*k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		
		free(W);

		PARDISO_numerical_factorization (A,A_pt,A_iparm,A_dparm,1);
	}
	else 
	{	  
		V  = mycalloc ("V of SPIKE_precond_setup'",k*k ,sizeof(double));
		W  = mycalloc ("W of SPIKE_precond_setup'",k*k ,sizeof(double));
		WT = mycalloc ("WT of SPIKE_precond_setup'",k*k ,sizeof(double));

		PARDISO_bottom_tips (A,&fB[0][0],V,k);
		PARDISO_top_tips    (A,&fC[0][0],W,k);

		MPI_Recv (WT,k*k,MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Isend(W ,k*k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);

		PARALLEL_REDUCED_system (V,WT,k,S,RANK, tag);
		
		free(V);
		free(W);
		free(WT);
		
		PARDISO_numerical_factorization (A,A_pt,A_iparm,A_dparm,1);
		PARDISO_numerical_factorization (S,S_pt,S_iparm,S_dparm,1);
	}
	
	SPIKE_precond (Parameters, MatrixData, FemStructs, &F[n_bef], &F[n_bef]);

	return 0;
}

/*----------------------------------------------------------------------------
 * Construct reduced system S for SPIKE preconditioner
 *--------------------------------------------------------------------------*/
void PARALLEL_REDUCED_system (double* V, double* W, int k, MAT* S, int RANK, int tag)
{
	int i,j,aa,ja;
	
	int nz = (2*k)+(2*k*k);
	int n  =  2*k;

	if (tag == 1){
		S->AA  = mycalloc("S->AA of PARALLEL_REDUCED_system'",nz ,sizeof(double));
		S->D   = mycalloc("S->D of PARALLEL_REDUCED_system'",n  ,sizeof(double));
		S->JA  = mycalloc("S->JA of PARALLEL_REDUCED_system'",nz ,sizeof(int   ));
		S->IA  = mycalloc("S->IA of PARALLEL_REDUCED_system'",n+1,sizeof(int   ));
		S->n   = n;
		S->m   = n;
		S->nz  = nz;
	}
	
	/* IA */
	for (i = 1; i < 2*k + 1; ++i)
		S->IA[i] = S->IA[i-1] + k + 1;

	/* AA, JA, g */
	aa = 0;
	ja = 0;
	for (i = 0; i < k; ++i)
	{
		S->AA[aa++] = 1.0;
		S->JA[ja++] = i;
		for (j = 0; j < k; ++j)
		{
			S->AA[aa++] = V[i*k + j];
			S->JA[ja++] = j + k;
		}
	}
	
	for (i = k; i < 2*k; ++i)
	{
		for (j = 0; j < k; ++j)
		{
			S->AA[aa++] = W[(i-k)*k + j];
			S->JA[ja++] = j;
		}
		S->AA[aa++] = 1.0;
		S->JA[ja++] = i;
	}
}


