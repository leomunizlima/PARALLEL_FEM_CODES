#include "../preconditioners.h"
#include "../pardiso.h"
#include <mpi.h>
#include "../../Allocation_Operations/allocations.h"

void MATRIX_matvec (MAT* A, double* u, double* p);

int SPIKE_precond(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *u, double *z)
{
	MAT *A		= MatrixData->A_SPIKE;
	MAT *B		= MatrixData->B_SPIKE;
	MAT *C		= MatrixData->C_SPIKE;
	MAT *S		= MatrixData->S_SPIKE;
	int n		= Parameters->neq;
	int k		= Parameters->k_SPIKE_block_size;
	void **A_pt 	= Parameters->PardisoVariables->A_pt;
	int *A_iparm 	= Parameters->PardisoVariables->A_iparm;
	double *A_dparm = Parameters->PardisoVariables->A_dparm;
	void **S_pt	= Parameters->PardisoVariables->S_pt; 
	int *S_iparm	= Parameters->PardisoVariables->S_iparm; 
	double *S_dparm	= Parameters->PardisoVariables->S_dparm; 
	double *f	= u; 
	double *y 	= z;
	int RANK 	= Parameters->RANK;
	int NPROC	= Parameters->NPROC;
	int i;
	MPI_Request request;
	
	double *g, *gT, *gS;
	double *x, *xB, *xT, *xTp;
	double *r1,*r2,*p;
	

	if (RANK == 0)
	{
		g  = mycalloc ("g of SPIKE_precond'",n  ,sizeof(double));
		gT = mycalloc ("gT of SPIKE_precond'",k  ,sizeof(double));
		gS = mycalloc ("gS of SPIKE_precond'",2*k,sizeof(double));
		x  = mycalloc ("x of SPIKE_precond'",2*k,sizeof(double));
		xB = mycalloc ("xB of SPIKE_precond'",k  ,sizeof(double));
		xT = mycalloc ("xT of SPIKE_precond'",k  ,sizeof(double));
		r1 = mycalloc ("r1 of SPIKE_precond'",n  ,sizeof(double));
		p  = mycalloc ("p of SPIKE_precond'",n  ,sizeof(double));
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,f,g,1);
		
		MPI_Recv(gT,k,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		memcpy(&gS[0],&g[n-k],k*sizeof(double));
		memcpy(&gS[k],&gT[0] ,k*sizeof(double));
		
		PARDISO_back_substitution (S,S_pt,S_iparm,S_dparm,gS,x,1);
		
		memcpy(&xB[0],&x[0],k*sizeof(double));
		memcpy(&xT[0],&x[k],k*sizeof(double));
		
		MPI_Isend(xB,k,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		
		MATRIX_matvec (B,xT,r1);
		for (i = 0; i < n; ++i)
			p[i] = f[i] - r1[i];
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,p,y,1);
	}
	else if (RANK == NPROC - 1)
	{
		g  = mycalloc ("g of SPIKE_precond'",n,sizeof(double));
		xB = mycalloc ("xB of SPIKE_precond'",k,sizeof(double));
		r1 = mycalloc ("r1 of SPIKE_precond'",n,sizeof(double));
		p  = mycalloc ("p of SPIKE_precond'",n,sizeof(double));
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,f,g,1);
		
		MPI_Isend(g, k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (xB,k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		MATRIX_matvec (C,xB,r1);
		for (i = 0; i < n; ++i)
			p[i] = f[i] - r1[i];
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,p,y,1);
	}
	else
	{
		g  = mycalloc ("xB of SPIKE_precond'",n  ,sizeof(double));
		gT = mycalloc ("gT of SPIKE_precond'",k  ,sizeof(double));
		gS = mycalloc ("gS of SPIKE_precond'",2*k,sizeof(double));
		x  = mycalloc ("x of SPIKE_precond'",2*k,sizeof(double));
		xB = mycalloc ("xB of SPIKE_precond'",k  ,sizeof(double));
		xT = mycalloc ("xT of SPIKE_precond'",k  ,sizeof(double));
		xTp = mycalloc("xTp of SPIKE_precond'",k  ,sizeof(double));
		r1 = mycalloc ("r1 of SPIKE_precond'",n  ,sizeof(double));
		r2 = mycalloc ("r2 of SPIKE_precond'",n  ,sizeof(double));
		p  = mycalloc ("p of SPIKE_precond'",n  ,sizeof(double));

		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,f,g,1);
		
		MPI_Isend(g, k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (gT,k,MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		memcpy(&gS[0],&g[n-k],k*sizeof(double));
		memcpy(&gS[k],&gT[0] ,k*sizeof(double));
		
		PARDISO_back_substitution (S,S_pt,S_iparm,S_dparm,gS,x,1);
		
		memcpy(&xB[0],&x[0],k*sizeof(double));
		memcpy(&xT[0],&x[k],k*sizeof(double));		
	
		MPI_Isend(xB ,k,MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (xTp,k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		MATRIX_matvec (C,xTp,r1);
		MATRIX_matvec (B,xT,r2);
		for (i = 0; i < n; ++i)
			p[i] = f[i] - r1[i] - r2[i];
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,p,y,1);
	}
	
//	MPI_Barrier(MPI_COMM_WORLD);
	
	if (RANK == 0)
	{
		free(g); 
		free(gT);
		free(gS);
		free(x);
		free(xB);
		free(xT);
		free(r1);
		free(p);
	}
	else if (RANK == NPROC - 1)
	{
		free(g);
		free(xB);
		free(r1);
		free(p);
	}
	else
	{
		free(g);  
		free(gT); 
		free(gS); 
		free(x);  
		free(xB); 
		free(xT); 
		free(xTp);
		free(r1); 
		free(r2); 
		free(p); 
	}

	

	return 0;
}

/*----------------------------------------------------------------------------
 * Compute the matrix-vector operation  p = Au
 *--------------------------------------------------------------------------*/
void MATRIX_matvec (MAT* A, double* u, double* p)
{
	int i, j, k1, k2;
	int n = A->n;
	double soma = 0;
	
	for (i = 0; i < n; ++i)
	{
		soma = 0;
		k1 = A->IA[i];
		k2 = A->IA[i+1]-1;
		for (j = k1; j <= k2; ++j)
			soma = soma + (A->AA[j]) * u[A->JA[j]];
		p[i] = soma;
	}  
}


