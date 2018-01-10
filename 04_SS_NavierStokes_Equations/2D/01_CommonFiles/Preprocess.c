#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

int Preprocess(int RANK, int NPROC, int narg, char **arguments, ParametersType **Parameters_out, MatrixDataType **MatrixData_out, FemStructsType **FemStructs_out, 
		FemFunctionsType **FemFunctions_out, FemOtherFunctionsType **FemOtherFunctions_out)
{
	int neq, nnodes, nel, I, J;
	int tag = 1; // Testing input error
	int neq_bef; //Number of equations on the previous partition 
	int NEQ; //Number of equations in subsequent partitions (RANK 0: d[2]; otherwise: d[RANK]-d[RANK-2])
	int nsend_bef, nsend_aft, nrecv_bef, nrecv_aft; //number of values to send and receive to RANK-1 and RANK+1
	int **lm, *lmaux;
	int size = NDOF*NNOEL;
	int size2 = size*size;
	double *F, *U, *f, *u, *Diag, *invDiag;
	char FileName[200], label[200];
	FILE *InFile;
	NodeType *Node;
	ElementType *Element;
	ParametersType *Parameters;
	MatrixDataType *MatrixData;
	FemStructsType *FemStructs;
	FemFunctionsType *FemFunctions;
	FemOtherFunctionsType *FemOtherFunctions;

	/* **************************************************************************************************************************** */
	//					Testing initial parameters
	/* **************************************************************************************************************************** */
	if (narg!=2)
	{
		if (RANK==0) printf("Use ./SSNavierStokesEEquations2D <Parameters file according to README>\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if (NPROC<2){
		printf("This code is developed for at least 2 MPI ranks\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	/* **************************************************************************************************************************** */


	/* **************************************************************************************************************************** */
	//			Reading parameters from problem setting file
	/* **************************************************************************************************************************** */
	Parameters = mycalloc("Parameters of 'Preprocess'", 1, sizeof(ParametersType));
	MatrixData   = mycalloc("MatrixData of 'Preprocess'",1,sizeof(MatrixDataType));
	FemStructs   = mycalloc("FemStructs of 'Preprocess'",1,sizeof(FemStructsType));
	FemFunctions = mycalloc("FemFunctions of 'Preprocess'",1,sizeof(FemFunctionsType));
	FemOtherFunctions = mycalloc("FemOtherFunctions of 'Preprocess'",1,sizeof(FemOtherFunctionsType));
	
	
	if (RANK == 0){
		InFile = myfopen(arguments[1], "r");
		tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->Experiments, label);
		tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->ProblemTitle, label);
		tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->ReynoldsNumber), label);
		tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->SolverTolerance), label);
		tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->LinearMaxIter), label);
		tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->NonLinearTolerance), label);
		tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->NonLinearMaxIter), label);
		tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->KrylovBasisVectorsQuantity), label);
		tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->Solver, label);
		tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->Preconditioner, label);
		tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->reordering, label);
		tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->MatrixVectorProductScheme, label);
		tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->StabilizationForm, label);
		tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->nNodes), label);
		tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->nEl), label);
		fclose(InFile);

		for (I=1; I<NPROC; I++) {
			MPI_Send(Parameters->Experiments, 			200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->ProblemTitle, 			200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->ReynoldsNumber), 		1, MPI_DOUBLE, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->SolverTolerance), 		1, MPI_DOUBLE, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->LinearMaxIter),	 		1,    MPI_INT, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->NonLinearTolerance), 		1, MPI_DOUBLE, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->NonLinearMaxIter),		1,    MPI_INT, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->KrylovBasisVectorsQuantity),	1,    MPI_INT, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->Solver,	 			200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->Preconditioner, 			200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->reordering, 			200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->MatrixVectorProductScheme,		200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->StabilizationForm,			200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->nNodes), 			1,    MPI_INT, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->nEl),	 			1,    MPI_INT, I, 0, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(Parameters->Experiments, 			200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->ProblemTitle, 			200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->ReynoldsNumber), 		1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->SolverTolerance), 		1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->LinearMaxIter),	 		1,    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->NonLinearTolerance), 		1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->NonLinearMaxIter),		1,    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->KrylovBasisVectorsQuantity),	1,    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->Solver,	 			200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->Preconditioner, 			200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->reordering, 			200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->MatrixVectorProductScheme,		200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->StabilizationForm,			200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->nNodes), 			1,    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->nEl),	 			1,    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}	
	/* **************************************************************************************************************************** */


	/* **************************************************************************************************************************** */
	//						Reading nodes
	/* **************************************************************************************************************************** */
	sprintf(FileName,"../02_mesh/%02d/RANK_%02d_%s_%d_%d.dat", NPROC, RANK, Parameters->ProblemTitle, Parameters->nNodes, Parameters->nEl);

	MPI_Barrier(MPI_COMM_WORLD);

	

	InFile = myfopen(FileName, "r");
	tag = fscanf(InFile, "%d", &nnodes);

	Node = (NodeType*) mycalloc("Node of 'Preprocess'", nnodes, sizeof(NodeType));
	for (I=0, neq = 0, nsend_bef = 0, nsend_aft = 0, nrecv_bef = 0, nrecv_aft = 0; I<nnodes; I++)
	{
		tag = fscanf(InFile, "%lf%lf", &(Node[I].x), &(Node[I].y));

		for (J=0; J<NDOF; J++)
		{
			tag = fscanf (InFile, "%d", &Node[I].Type[J]);

			if (Node[I].Type[J] == 0){
				Node[I].invP_id[J] = -1;
			}
			else {
				tag =  fscanf(InFile, "%d%d", &(Node[I].Send[J]), &(Node[I].invP_id[J]));
				if(Node[I].Type[J] == 1){
					neq++; 
					if (Node[I].Send[J] == 2 || Node[I].Send[J] == 23)
						nsend_bef++;
					if (Node[I].Send[J] == 3 || Node[I].Send[J] == 23)
						nsend_aft++;
				}else if (Node[I].Type[J] == 2)
					nrecv_bef++;
				else
					nrecv_aft++;
			}
		}
	}

	printf("nsend_bef=%d nsend_aft=%d nrecv_bef=%d nrecv_aft=%d (RANK=%d)\n",nsend_bef,nsend_aft,nrecv_bef,nrecv_aft,RANK);

	/* **************************************************************************************************************************** */


	/* **************************************************************************************************************************** */
	//           				Reading connection mesh
	/* **************************************************************************************************************************** */
	tag = fscanf(InFile, "%d", &nel);
	Element = (ElementType*) mycalloc("Element of 'Preprocess'", nel, sizeof(ElementType));
	for (I = 0; I < nel; I++)
		tag = fscanf(InFile, "%d%d%d%d", &(Element[I].Vertex[0]), &(Element[I].Vertex[1]), &(Element[I].Vertex[2]), &(Element[I].Type));
	fclose(InFile);
	
	/* **************************************************************************************************************************** */

	/**********************************************************************************************/
	//			Setting buffers to send and receive data to RANK-1 and RANK+1
	/*********************************************************************************************/
	
	if (RANK != 0){
		int *IdSend_bef, *IdRecv_bef, *IdRecvAux_bef, *IdDegreeRecvAux_bef;

		FemStructs->SendBuffer_bef[0] = mycalloc("SendBuffer_bef of 'Preprocess'",nsend_bef,sizeof(double));	
		FemStructs->SendBuffer_bef[1] = mycalloc("SendBuffer_bef of 'Preprocess'",nsend_bef,sizeof(double));	
		FemStructs->RecvBuffer_bef[0] = mycalloc("RecvBuffer_bef of 'Preprocess'",nrecv_bef,sizeof(double));	
		FemStructs->RecvBuffer_bef[1] = mycalloc("RecvBuffer_bef of 'Preprocess'",nrecv_bef,sizeof(double));	
		IdSend_bef          = mycalloc("IdSend_bef of 'Preprocess'",nsend_bef,sizeof(int));
		IdRecv_bef          = mycalloc("IdRecv_bef of 'Preprocess'",nrecv_bef,sizeof(int));
		IdRecvAux_bef       = mycalloc("IdRecvAux_bef of 'Preprocess'",nrecv_bef,sizeof(int));
		IdDegreeRecvAux_bef = mycalloc("IdRecvAux_bef of 'Preprocess'",nrecv_bef,sizeof(int));
		ARRAYType *SortIdRecv_bef  = mycalloc("SortIdRecv_bef of 'Preprocess'",nrecv_bef,sizeof(ARRAYType));

		FemStructs->IdSend_bef = IdSend_bef;
		FemStructs->IdRecv_bef = IdRecv_bef;
	
		int i_bef, j_bef;
			
		for (I=0, i_bef=0, j_bef=0; I<nnodes; I++){
			for (J=0; J<NDOF; J++){
				if (Node[I].Send[J] == 2 || Node[I].Send[J] == 23)
					IdSend_bef[i_bef++] = Node[I].invP_id[J];
				if (Node[I].Type[J] == 2){
					IdRecv_bef[j_bef] = Node[I].invP_id[J];
					IdRecvAux_bef[j_bef] = I;
					IdDegreeRecvAux_bef[j_bef] = J;
					j_bef++;
				}
			}
		}
		//-----------------------------------------------------------------------------------------//
		//  		Adjustment node indexes that came from RANK-1
		//-----------------------------------------------------------------------------------------//

		for (I=0; I<nrecv_bef; I++){
			SortIdRecv_bef[I].array1 = IdRecv_bef[I];
			SortIdRecv_bef[I].array2 = I;
		}

		qsort(SortIdRecv_bef,nrecv_bef,sizeof(ARRAYType),COMPARE_array_in_Preprocess);
	
		for (I=0; I<nrecv_bef; I++){
			while (SortIdRecv_bef[I].array1 + nrecv_bef < I) 
				SortIdRecv_bef[I].array1++;
		}

		for (I=0; I<nrecv_bef; I++){
			IdRecv_bef[SortIdRecv_bef[I].array2] = SortIdRecv_bef[I].array1;
			Node[IdRecvAux_bef[SortIdRecv_bef[I].array2]].invP_id[IdDegreeRecvAux_bef[SortIdRecv_bef[I].array2]] = SortIdRecv_bef[I].array1;
		}

		free(SortIdRecv_bef);
		free(IdRecvAux_bef);
		free(IdDegreeRecvAux_bef);

	}

	if (RANK != NPROC-1){

		int *IdSend_aft, *IdRecv_aft, *IdRecvAux_aft, *IdDegreeRecvAux_aft;

		FemStructs->SendBuffer_aft[0] = mycalloc("SendBuffer_aft of 'Preprocess'",nsend_aft,sizeof(double));	
		FemStructs->SendBuffer_aft[1] = mycalloc("SendBuffer_aft of 'Preprocess'",nsend_aft,sizeof(double));	
		FemStructs->RecvBuffer_aft[0] = mycalloc("RecvBuffer_aft of 'Preprocess'",nrecv_aft,sizeof(double));	
		FemStructs->RecvBuffer_aft[1] = mycalloc("RecvBuffer_aft of 'Preprocess'",nrecv_aft,sizeof(double));	
		IdSend_aft                 = mycalloc("IdSend_aft of 'Preprocess'",nsend_aft,sizeof(int));
		IdRecv_aft                 = mycalloc("IdRecv_aft of 'Preprocess'",nrecv_aft,sizeof(int));
		IdRecvAux_aft              = mycalloc("IdRecvAux_aft of 'Preprocess'",nrecv_aft,sizeof(int));
		IdDegreeRecvAux_aft        = mycalloc("IdRecvAux_aft of 'Preprocess'",nrecv_aft,sizeof(int));
		ARRAYType *SortIdRecv_aft  = mycalloc("SortIdRecv_aft of 'Preprocess'",nrecv_aft,sizeof(ARRAYType));

		FemStructs->IdSend_aft = IdSend_aft;
		FemStructs->IdRecv_aft = IdRecv_aft;
	

		int i_aft, j_aft;

		for (I=0, i_aft=0, j_aft=0; I<nnodes; I++){
			for (J=0; J<NDOF; J++){
				if (Node[I].Send[J] == 3 || Node[I].Send[J] == 23)
					IdSend_aft[i_aft++] = Node[I].invP_id[J];
				if (Node[I].Type[J] == 3){
					IdRecv_aft[j_aft] = Node[I].invP_id[J];
					IdRecvAux_aft[j_aft] = I;
					IdDegreeRecvAux_aft[j_aft] = J;
					j_aft++;
				}
			}
		}

		//-----------------------------------------------------------------------------------------//
		//  		Adjustment node indexes that came from RANK+1
		//-----------------------------------------------------------------------------------------//

		for (I=0; I<nrecv_aft; I++){
			SortIdRecv_aft[I].array1 = IdRecv_aft[I];
			SortIdRecv_aft[I].array2 = I;
		}

		qsort(SortIdRecv_aft,nrecv_aft,sizeof(ARRAYType),COMPARE_array_in_Preprocess);
	
		for (I=0; I<nrecv_aft; I++){
			while (SortIdRecv_aft[I].array1 > I + neq) 
				SortIdRecv_aft[I].array1--;
		}

		for (I=0; I<nrecv_aft; I++){
			IdRecv_aft[SortIdRecv_aft[I].array2] = SortIdRecv_aft[I].array1;
			Node[IdRecvAux_aft[SortIdRecv_aft[I].array2]].invP_id[IdDegreeRecvAux_aft[SortIdRecv_aft[I].array2]] = SortIdRecv_aft[I].array1;
		}

		free(SortIdRecv_aft);
		free(IdRecvAux_aft);
		free(IdDegreeRecvAux_aft);


	}

	/* **************************************************************************************************************************** */
	//			          Memory allocations and Store strategies 
	/* **************************************************************************************************************************** */

	// Some variable initializations
	
	NEQ = nrecv_bef + neq + nrecv_aft;
	neq_bef = nrecv_bef;
	
	printf("neq_bef=%d neq=%d NEQ=%d nel=%d (RANK: %d)\n",neq_bef,neq,NEQ,nel,RANK);

	Parameters->neq = neq;
	Parameters->neq_bef = neq_bef;
	Parameters->nnodes = nnodes;
	Parameters->nel = nel;
	Parameters->NEQ = NEQ;
	Parameters->RANK = RANK;
	Parameters->NPROC = NPROC;
	Parameters->nsend_bef = nsend_bef;
	Parameters->nsend_aft = nsend_aft;
	Parameters->nrecv_bef = nrecv_bef;
	Parameters->nrecv_aft = nrecv_aft;

	F = (double *) mycalloc("F of 'Preprocess'", NEQ+1, sizeof(double));  
	U = (double *) mycalloc("U of 'Preprocess'", NEQ+1, sizeof(double)); // U has size NEQ+1: u_(RANK-1) + u_RANK + u_(RANK+1) + 1
	u = &U[neq_bef];
	f = &F[neq_bef];
	Diag = (double*) mycalloc("Diag of 'Preprocess'",NEQ+1, sizeof(double));
	invDiag = (double*) mycalloc("invDiag of 'Preprocess'",NEQ+1, sizeof(double));

	lm = (int**) mycalloc("lm of 'Preprocess'", nel, sizeof(int*));
	lmaux = (int*) mycalloc("lmaux of 'Preprocess'", nel*size, sizeof(int));
	for (I = 0; I < nel; I++)
		lm[I] = &lmaux[I*size];
		
	
	//Configuring equation according to variables and boundary conditions
	Fill_LM(Parameters, lm, Node, Element);
	
	FemStructs->Node = Node;
	FemStructs->Element = Element;
	FemStructs->lm = lm;
	FemStructs->lmaux = lmaux;
	FemStructs->F = F;
	FemStructs->U = U;
	FemStructs->f = f;
	FemStructs->u = u;

	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE") == 0){		
		double **A, *Aaux;

		A = (double**) mycalloc("K of 'Preprocess'", nel, sizeof(double*));
		Aaux = (double*) mycalloc("Kaux of 'Preprocess'", nel*size2,sizeof(double));

		for (I = 0; I < nel; I++){
			A[I] = &Aaux[I*size2];
		}
		
		MatrixData->A = A;
		MatrixData->Aaux = Aaux;
		
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") == 0){
		double *AA, *AA_bef, *AA_aft;

		csr_Initialization(Parameters, MatrixData, FemStructs);
		
		AA = (double*) mycalloc("AA of 'Preprocess'", Parameters->nnzero+1, sizeof(double));
		AA_bef = (double*) mycalloc("AA of 'Preprocess'", Parameters->nnzero_bef+1, sizeof(double));
		AA_aft = (double*) mycalloc("AA of 'Preprocess'", Parameters->nnzero_aft+1, sizeof(double));
		
		MatrixData->AA = AA;
		MatrixData->AA_bef = AA_bef;
		MatrixData->AA_aft = AA_aft;
	}
	else{
		if (RANK ==0) printf("Matrix vector product scheme is not defined correctly!\n\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	/* **************************************************************************************************************************** */

	MatrixData->Diag = Diag;
	MatrixData->invDiag = invDiag;

	*Parameters_out = Parameters;
	*MatrixData_out = MatrixData;
	*FemStructs_out = FemStructs;
	*FemFunctions_out = FemFunctions;
	*FemOtherFunctions_out = FemOtherFunctions;

	if (tag<0){
		if (RANK == 0) printf ("Error in some parameter\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	
	return 0;

}


