#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

int Preprocess(int RANK, int NPROC,  int narg, char **arguments, ParametersType **Parameters_out, MatrixDataType **MatrixData_out, FemStructsType **FemStructs_out, 
		FemFunctionsType **FemFunctions_out, FemOtherFunctionsType **FemOtherFunctions_out) 
{
	int neq; //Number of equations on the partition 
	int neq_bef; //Number of equations on the previous partition 
	int NEQ; //Number of equations in subsequent partitions (RANK 0: d[2]; otherwise: d[RANK]-d[RANK-2])
	int nsend_bef, nsend_aft, nrecv_bef, nrecv_aft; //number of values to send and receive to RANK-1 and RANK+1
	int I, nel, nnodes;
	int **lm, *lmaux;
	int *Perm, *invPerm;
	int size = NDOF*NNOEL;
	int size2 = size*size;
	double *F, *f, *U, *u, *Diag, *invDiag;
	NodeType *Node;
	ElementType *Element;
	char FileName[300];
	char label[300];
	int tag = 1; // Testing input error
	FILE *InFile;
	ParametersType *Parameters;
	MatrixDataType *MatrixData;
	FemStructsType *FemStructs;
	FemFunctionsType *FemFunctions;
	FemOtherFunctionsType *FemOtherFunctions;


	//-----------------------------------Example with 3 MPI ranks------------------------------------
	//		| 				|           |	 |	--> Ai      means a matrix square block with size neq_i 
	//		|   A11	     A1_aft	    	|	    | F1 | 	--> Ai_bef  means a matrix block with size neq_i x neq_{i-1}
	//		|				|	    |	 |	--> Ai_aft  means a matrix block with size neq_i x neq_{i+1}
        //	A =	| A2_bef      A2        A2_aft 	| and	F = | F2 |	--> NEQ     means size neq_{i-1} + neq_i + neq_{i+1}
	//		|				|	    |	 |
	//		|            A3_bef	  A3   	|	    | F3 |
        //		|				|	    |	 |
        //
	//--------------------------------------------------------------------------------

	/*****************************************************************************************************************************************/
	//							Testing initial parameters
	/******************************************************************************************************************************************/
	if (narg!=4)
	{
		if (RANK==0) printf("Use ./SSTranspEquation <Parameters file according to README> <path to input files> <path to output files>\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if (NPROC<2){
		printf("This code is developed for at least 2 MPI ranks\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	/***************************************************************************/
//			Reading parameters from problem setting file
	/***************************************************************************/
	Parameters   = mycalloc("Parameters of 'Preprocess'",1,sizeof(ParametersType));	
	MatrixData   = mycalloc("MatrixData of 'Preprocess'",1,sizeof(MatrixDataType));	
	FemStructs   = mycalloc("FemStructs of 'Preprocess'",1,sizeof(FemStructsType));
	FemFunctions = mycalloc("FemFunctions of 'Preprocess'",1,sizeof(FemFunctionsType));
	FemOtherFunctions = mycalloc("FemOtherFunctions of 'Preprocess'",1,sizeof(FemOtherFunctionsType));
	Parameters->arguments = arguments;
	
	if (RANK==0){
		InFile = myfopen(arguments[1], "r");
		tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->ProblemTitle, label);
		tag = fscanf(InFile, "%lf\t:%[^\n]\n",&(Parameters->SolverTolerance), label);
		tag = fscanf(InFile, "%lf\t:%[^\n]\n",&(Parameters->NonLinearTolerance), label);
		tag = fscanf(InFile, "%d\t:%[^\n]\n", &(Parameters->LinearMaxIter), label);
		tag = fscanf(InFile, "%d\t:%[^\n]\n", &(Parameters->KrylovBasisVectorsQuantity), label);
		tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->Solver, label);
		tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->Preconditioner, label);
		tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->Scaling, label);
		tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->reordering, label);
		tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->MatrixVectorProductScheme, label);
		tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->StabilizationForm, label);
		tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->ShockCapture, label);
		tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->h_Shock, label);
		tag = fscanf(InFile, "%d\t:%[^\n]\n", &(Parameters->nNodes), label);
		tag = fscanf(InFile, "%d\t:%[^\n]\n", &(Parameters->nEl), label);
		fclose(InFile);
	

		for (I=1; I<NPROC;I++){
			MPI_Send(Parameters->ProblemTitle,                  200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->SolverTolerance),            1, MPI_DOUBLE, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->NonLinearTolerance),         1, MPI_DOUBLE, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->LinearMaxIter),              1,    MPI_INT, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->KrylovBasisVectorsQuantity), 1,    MPI_INT, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->Solver,                        200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->Preconditioner,                200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->Scaling,                200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->reordering,                    200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->MatrixVectorProductScheme,     200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->StabilizationForm,             200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->ShockCapture,                  200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(Parameters->h_Shock,                       200, MPI_CHAR, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->nNodes),                     1,    MPI_INT, I, 0, MPI_COMM_WORLD);
			MPI_Send(&(Parameters->nEl),                        1,    MPI_INT, I, 0, MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Recv(Parameters->ProblemTitle,                  200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->SolverTolerance),            1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->NonLinearTolerance),         1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->LinearMaxIter),              1,    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->KrylovBasisVectorsQuantity), 1,    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->Solver,                        200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->Preconditioner,                200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->Scaling,                200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->reordering,                    200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->MatrixVectorProductScheme,     200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->StabilizationForm,             200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->ShockCapture,                  200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Parameters->h_Shock,                       200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->nNodes),                     1,    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(Parameters->nEl),                        1,    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	/*****************************************************************************/


	/*****************************************************************************/
	//				Reading nodes	
	/*****************************************************************************/		
	sprintf(FileName,"%s/%02d/RANK_%02d_%s_%d_%d.dat", arguments[2], NPROC, RANK, Parameters->ProblemTitle, Parameters->nNodes, Parameters->nEl);
	InFile = myfopen(FileName, "r");
	tag =  fscanf(InFile, "%d", &nnodes);


	Node = (NodeType*) mycalloc("Node of 'Preprocess'",nnodes, sizeof(NodeType));
	for (I=0, neq = 0, nsend_bef = 0, nsend_aft = 0, nrecv_bef = 0, nrecv_aft = 0; I<nnodes; I++)
	{
		tag = fscanf(InFile, "%lf%lf%d", &(Node[I].x), &(Node[I].y), &(Node[I].Type));

		if (Node[I].Type == 0){
			Node[I].invP_id = -1;
		}
		else {
			tag =  fscanf(InFile, "%d%d", &(Node[I].Send), &(Node[I].invP_id));
			if(Node[I].Type == 1){
				neq++; 
				if (Node[I].Send == 2 || Node[I].Send == 23)
					nsend_bef++;
				if (Node[I].Send == 3 || Node[I].Send == 23)
					nsend_aft++;
			}else if (Node[I].Type == 2)
				nrecv_bef++;
			else
				nrecv_aft++;
		}
	}

	printf("nsend_bef=%d nsend_aft=%d nrecv_bef=%d nrecv_aft=%d (RANK=%d)\n",nsend_bef,nsend_aft,nrecv_bef,nrecv_aft,RANK);

	/*****************************************************************************/


	/************************************************************************************************************/	
	//						Reading connection mesh
	/************************************************************************************************************/	
	tag = fscanf(InFile, "%d", &nel);
	Element = (ElementType*) mycalloc("Element of 'Preprocess'",nel, sizeof(ElementType));
	for (I=0; I<nel; I++)
		tag = fscanf(InFile, "%d%d%d", &(Element[I].Vertex[0]),  &(Element[I].Vertex[1]), &(Element[I].Vertex[2]));
	
	fclose(InFile);
	/*************************************************************************************************************/
	
	/**********************************************************************************************/
	//			Setting buffers to send and receive data to RANK-1 and RANK+1
	/*********************************************************************************************/
	int *IdSend_bef, *IdSend_aft, *IdRecv_bef, *IdRecv_aft, *IdRecvAux_bef, *IdRecvAux_aft;
	
	if (RANK != 0){
		FemStructs->SendBuffer_bef = mycalloc("SendBuffer_bef of 'Preprocess'",nsend_bef,sizeof(double));	
		FemStructs->RecvBuffer_bef = mycalloc("RecvBuffer_bef of 'Preprocess'",nrecv_bef,sizeof(double));	

		IdSend_bef          = mycalloc("IdSend_bef of 'Preprocess'",nsend_bef,sizeof(int));
		IdRecv_bef          = mycalloc("IdRecv_bef of 'Preprocess'",nrecv_bef,sizeof(int));
		IdRecvAux_bef       = mycalloc("IdRecvAux_bef of 'Preprocess'",nrecv_bef,sizeof(int));
		ARRAYType *SortIdRecv_bef  = mycalloc("SortIdRecv_bef of 'Preprocess'",nrecv_bef,sizeof(ARRAYType));

		FemStructs->IdSend_bef = IdSend_bef;
		FemStructs->IdRecv_bef = IdRecv_bef;
	
		int i_bef, j_bef;
			
		for (I=0, i_bef=0, j_bef=0; I<nnodes; I++){
			if (Node[I].Send == 2 || Node[I].Send == 23)
				IdSend_bef[i_bef++] = Node[I].invP_id;
			if (Node[I].Type == 2){
				IdRecv_bef[j_bef] = Node[I].invP_id;
				IdRecvAux_bef[j_bef] = I;
				j_bef++;
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
			Node[IdRecvAux_bef[SortIdRecv_bef[I].array2]].invP_id = SortIdRecv_bef[I].array1;
		}

		free(SortIdRecv_bef);
		free(IdRecvAux_bef);

	}

	if (RANK != NPROC-1){
		FemStructs->SendBuffer_aft = mycalloc("SendBuffer_aft of 'Preprocess'",nsend_aft,sizeof(double));	
		FemStructs->RecvBuffer_aft = mycalloc("RecvBuffer_aft of 'Preprocess'",nrecv_aft,sizeof(double));	
		IdSend_aft                 = mycalloc("IdSend_aft of 'Preprocess'",nsend_aft,sizeof(int));
		IdRecv_aft                 = mycalloc("IdRecv_aft of 'Preprocess'",nrecv_aft,sizeof(int));
		IdRecvAux_aft              = mycalloc("IdRecvAux_aft of 'Preprocess'",nrecv_aft,sizeof(int));
		ARRAYType *SortIdRecv_aft  = mycalloc("SortIdRecv_aft of 'Preprocess'",nrecv_aft,sizeof(ARRAYType));

		FemStructs->IdSend_aft = IdSend_aft;
		FemStructs->IdRecv_aft = IdRecv_aft;
	

		int i_aft, j_aft;

		for (I=0, i_aft=0, j_aft=0; I<nnodes; I++){
			if (Node[I].Send == 3 || Node[I].Send == 23)
				IdSend_aft[i_aft++] = Node[I].invP_id;
			if (Node[I].Type == 3){
				IdRecv_aft[j_aft] = Node[I].invP_id;
				IdRecvAux_aft[j_aft] = I;
				j_aft++;
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
			Node[IdRecvAux_aft[SortIdRecv_aft[I].array2]].invP_id = SortIdRecv_aft[I].array1;
		}

		free(SortIdRecv_aft);
		free(IdRecvAux_aft);


	}

	/*****************************************************************************/
	//                      Memory allocations and Store strategies 
        /*****************************************************************************/

	// Some variable inicializations

	NEQ = nrecv_bef + neq + nrecv_aft;
	neq_bef = nrecv_bef;
	
	printf("neq_bef=%d neq=%d NEQ=%d nel=%d (RANK: %d)\n",neq_bef,neq,NEQ,nel,RANK);

	Parameters->neq = neq;
	Parameters->neq_bef = neq_bef;
	Parameters->nel = nel;
	Parameters->nnodes = nnodes;
	Parameters->NEQ = NEQ;
	Parameters->RANK = RANK;
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

	lm = (int**) mycalloc("lm of 'Preprocess'",nel, sizeof(int*));
	lmaux = (int*) mycalloc("lmaux of 'Preprocess'",nel*size, sizeof(int));
	for (I=0;I<nel;I++)
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
	
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){		
	
		//Reducing bandwidth to decrease 'cache miss' inside matrix vector product per MPI rank
		Perm    = (int *) mycalloc("Perm of 'Preprocess'",NEQ+1, sizeof(int)); 		
		invPerm    = (int *) mycalloc("Perm of 'Preprocess'",NEQ+1, sizeof(int)); 		
		
		Permutation_of_LM(Parameters, FemStructs, Perm, invPerm);

		double **K, *Kaux;		

		K = (double**) mycalloc("K of 'Preprocess'", nel, sizeof(double*));
		Kaux = (double*) mycalloc("Kaux of 'Preprocess'", nel*size2,sizeof(double));

		for (I=0;I<nel;I++){
			K[I] = &Kaux[I*size2];
		}
		MatrixData->A = K;
		MatrixData->Aaux = Kaux;
		MatrixData->Perm = Perm;
		MatrixData->invPerm = invPerm;
	} 
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){		
	
		double *K, *K_bef, *K_aft;

		csr_Initialization(Parameters, MatrixData, FemStructs);

		K = (double*) mycalloc("K of 'Preprocess'", Parameters->nnzero+1, sizeof(double));
		K_bef = (double*) mycalloc("K of 'Preprocess'", Parameters->nnzero_bef+1, sizeof(double));
		K_aft = (double*) mycalloc("K of 'Preprocess'", Parameters->nnzero_aft+1, sizeof(double));

		MatrixData->AA = K;
		MatrixData->AA_bef = K_bef;
		MatrixData->AA_aft = K_aft;
	} 
	else{
		if (RANK==0) printf("Matrix vector product scheme not defined!\n\n");
		MPI_Finalize();
	}
	/************************************************************************************/

	MatrixData->Diag = Diag;
        MatrixData->invDiag = invDiag;
	Parameters->RANK = RANK;	
	Parameters->NPROC = NPROC;
	*Parameters_out = Parameters;
	*MatrixData_out = MatrixData;
	*FemStructs_out = FemStructs;
	*FemFunctions_out = FemFunctions;
	*FemOtherFunctions_out = FemOtherFunctions;

	if (tag<0){
		 if (RANK==0) printf ("Error in some parameter\n");
		 MPI_Finalize();
	}

	return 0;

}


