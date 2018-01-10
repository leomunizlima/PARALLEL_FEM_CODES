#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

void SPIKEMAT_clean(MAT *A);

int Postprocess(int RANK, int NPROC, ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int nEq;
	char FileName[200];
	FILE *OutFile;

	MPI_Allreduce(&Parameters->neq, &nEq, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	/*************************************************************/
	//		Paraview output to file
	/************************************************************/	
	Paraview_Output(RANK, NPROC, Parameters, FemStructs, FemFunctions);
	/*************************************************************/


	/****************************************************************************************/
		// 			Printing final result		
	/****************************************************************************************/
	if (RANK == 0){
		sprintf(FileName, "%s/%02d/%s_%s_%s_%s_%s_%s_N%d_E%d.dat",Parameters->arguments[3],NPROC,Parameters->ProblemTitle,Parameters->StabilizationForm,
			Parameters->ShockCapture,Parameters->h_Shock, Parameters->MatrixVectorProductScheme,Parameters->Preconditioner, Parameters->nNodes,Parameters->nEl); 	
		OutFile = myfopen(FileName,"w");
		fprintf(OutFile, "\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
		fprintf(OutFile, "Problem Title: Steady State %s\n", Parameters->ProblemTitle);
		fprintf(OutFile, "Number of nodes: %d\n", Parameters->nNodes);
		fprintf(OutFile, "Number of elements: %d\n", Parameters->nEl);
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
			fprintf(OutFile, "Number of edges: %d\n", Parameters->nedge);
		}	
		fprintf(OutFile, "Number of equations: %d\n", nEq);
		fprintf(OutFile, "Stabilization form used: %s\n", Parameters->StabilizationForm);
		fprintf(OutFile, "Shock capture used: %s\n", Parameters->ShockCapture);
		fprintf(OutFile, "Parameter h for shock capture used: %s\n", Parameters->h_Shock);
		fprintf(OutFile, "Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
			fprintf(OutFile,"Reordering: %s (bandwidth before: %d) (bandwidth after: %d)\n", Parameters->reordering, 
			Parameters->bandwidth_bef, Parameters->bandwidth_aft);
		}
		fprintf(OutFile, "Solver used: %s\n", Parameters->Solver);
		fprintf(OutFile, "Preconditioner used: %s\n", Parameters->Preconditioner);
		fprintf(OutFile, "Scaling used: %s\n", Parameters->Scaling);
		fprintf(OutFile, "Solver tolerance used: %E\n", Parameters->SolverTolerance);
		fprintf(OutFile, "Non linear tolerance used: %E\n", Parameters->NonLinearTolerance);
		fprintf(OutFile, "Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);		
		fprintf(OutFile, "Number of MPI ranks: %d\n", NPROC);
		fprintf(OutFile, "\n========================================================================\n\n");
		fclose(OutFile);
        
		printf("\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
		printf("Problem Title: Steady State %s\n", Parameters->ProblemTitle);
		printf("Number of nodes: %d\n", Parameters->nNodes);
		printf("Number of elements: %d\n", Parameters->nEl);
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
			printf("Number of edges: %d\n", Parameters->nedge);
		}	
		printf("Number of equations: %d\n", nEq);
		printf("Stabilization form used: %s\n", Parameters->StabilizationForm);
		printf("Shock capture used: %s\n", Parameters->ShockCapture);
		printf("Parameter h for shock capture used: %s\n", Parameters->h_Shock);
		printf("Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
		printf("Solver used: %s\n", Parameters->Solver);
		printf("Preconditioner used: %s\n", Parameters->Preconditioner);
		printf("Scaling used: %s\n", Parameters->Scaling);
		printf("Solver tolerance used: %E\n", Parameters->SolverTolerance);
		printf("Non linear tolerance used: %E\n", Parameters->NonLinearTolerance);
		printf("Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);		
		printf("Number of MPI ranks: %d\n", NPROC);
		printf("\n========================================================================\n\n");
        
		/****************************************************************************************/
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		printf("Reordering: %s (bandwidth before: %d) (bandwidth after: %d) RANK: %d \n", Parameters->reordering, 
		Parameters->bandwidth_bef, Parameters->bandwidth_aft, RANK);
	}

	/***************************************************************************************/
	//				Memory deallocation	
	/**************************************************************************************/
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){		
		free(MatrixData->A);
		free(MatrixData->Aaux);
		free(MatrixData->Perm);
		free(MatrixData->invPerm);
		free(FemStructs->lm);
		free(FemStructs->lmaux);
	}

/*	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){		
		free(MatrixData->A[0]);
		free(MatrixData->Aaux[0]);
		for (I = 0; I < nel; I++){
			free(MatrixData->Scheme_by_Element[I]);
			free(MatrixData->order[I]);
		}
		free(MatrixData->Scheme_by_Element);
		free(MatrixData->order);
		free(lm);
		free(lmaux);

	}
*/
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){		
		free(MatrixData->AA);
		free(MatrixData->IA);
		free(MatrixData->JA);
		free(MatrixData->Perm);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
		if (Parameters->neq_bef>0){
			free(MatrixData->AA_bef);
			free(MatrixData->IA_bef);
			free(MatrixData->JA_bef);
			free(MatrixData->aux_bef);
		}
		if (Parameters->NEQ-Parameters->neq-Parameters->neq_bef>0){
			free(MatrixData->AA_aft);
			free(MatrixData->IA_aft);
			free(MatrixData->JA_aft);
			free(MatrixData->aux_aft);
		}
		int I;
		for (I = 0; I < Parameters->nel; I++){
			free(MatrixData->Scheme_by_Element[I]);
			free(MatrixData->Scheme_by_Element_bef[I]);
			free(MatrixData->Scheme_by_Element_aft[I]);
		}
		free(MatrixData->Scheme_by_Element);
		free(MatrixData->Scheme_by_Element_bef);
		free(MatrixData->Scheme_by_Element_aft);
		if (strncmp(Parameters->Preconditioner,"ILU",3)==0){
			SPARILU_clean(MatrixData->ILUp);
			SPARMAT_clean(MatrixData->mat);
			free(MatrixData->Ailu);
		}
		if (strncmp(Parameters->Preconditioner,"SPIKE",5)==0){
			if (RANK != 0){
				free(MatrixData->fullCaux);
				free(MatrixData->fullC);
				SPIKEMAT_clean (MatrixData->C_SPIKE);
			}

			if (RANK != NPROC -1){
				free(MatrixData->fullBaux);
				free(MatrixData->fullB);
				SPIKEMAT_clean (MatrixData->B_SPIKE);
				SPIKEMAT_clean (MatrixData->S_SPIKE);
			}
			free(MatrixData->A_SPIKE);
			free(Parameters->PardisoVariables);
		}
	}

	free(MatrixData->Diag);
	free(MatrixData->invDiag);
	free(MatrixData);
	free(FemStructs->Node);
	free(FemStructs->Element);
	free(FemStructs->F);
	free(FemStructs->U);
	if (RANK != 0){
		free(FemStructs->SendBuffer_bef);
		free(FemStructs->RecvBuffer_bef);
		free(FemStructs->IdSend_bef);    	
		free(FemStructs->IdRecv_bef);
	}
	if (RANK != NPROC-1){
		free(FemStructs->SendBuffer_aft);
		free(FemStructs->RecvBuffer_aft);
		free(FemStructs->IdSend_aft);    	
		free(FemStructs->IdRecv_aft);
	}    
	
	free(FemStructs);
	free(FemFunctions);
	free(FemOtherFunctions);

	/***************************************************************************************/	
	

	
	return 0;
}

void SPIKEMAT_clean(MAT *A)
{
	free(A->AA);
	free(A->D);
	free(A->JA);
	free(A->IA);
	free(A);
}


