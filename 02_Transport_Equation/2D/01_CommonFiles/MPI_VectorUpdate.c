#include "TranspEquation.h"

void MPI_VectorUpdate(int RANK, int NPROC, ParametersType *Parameters, FemStructsType *FemStructs, double *P, int tag)
{
	MPI_Request request;
	int I;
	int ns_bef = Parameters->nsend_bef;
	int ns_aft = Parameters->nsend_aft;
	int nr_bef = Parameters->nrecv_bef;
	int nr_aft = Parameters->nrecv_aft;
	double *SBuff_bef = FemStructs->SendBuffer_bef[tag]; 
	double *SBuff_aft = FemStructs->SendBuffer_aft[tag];
	double *RBuff_bef = FemStructs->RecvBuffer_bef[tag];
	double *RBuff_aft = FemStructs->RecvBuffer_aft[tag];
	int *IdS_bef = FemStructs->IdSend_bef; 
	int *IdS_aft = FemStructs->IdSend_aft; 
	int *IdR_bef = FemStructs->IdRecv_bef; 
	int *IdR_aft = FemStructs->IdRecv_aft;

	for (I=0; I<ns_bef; I++)
		SBuff_bef[I] = P[IdS_bef[I]];
	for (I=0; I<ns_aft; I++)
		SBuff_aft[I] = P[IdS_aft[I]];

	if (RANK == 0)
	{	  
		MPI_Isend(SBuff_aft,ns_aft,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (RBuff_aft,nr_aft,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	else if (RANK == NPROC - 1)
	{
		MPI_Isend(SBuff_bef,ns_bef,MPI_DOUBLE,NPROC-2,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (RBuff_bef,nr_bef,MPI_DOUBLE,NPROC-2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}

	else
	{
		MPI_Isend(SBuff_bef,ns_bef,MPI_DOUBLE,RANK-1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Isend(SBuff_aft,ns_aft,MPI_DOUBLE,RANK+1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (RBuff_bef,nr_bef,MPI_DOUBLE,RANK-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv (RBuff_aft,nr_aft,MPI_DOUBLE,RANK+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}

	for (I=0; I<nr_bef; I++)
		P[IdR_bef[I]] = RBuff_bef[I];
	for (I=0; I<nr_aft; I++)
		P[IdR_aft[I]] = RBuff_aft[I]; 


}


