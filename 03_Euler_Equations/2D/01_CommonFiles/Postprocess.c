#include "EulerEquations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

void SPIKEMAT_clean(MAT *A);

int Postprocess(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int nEq;
	char FileName[200];
	FILE *OutFile;
	int RANK,  NPROC;

	RANK = Parameters->RANK;
	NPROC = Parameters->NPROC;

	MPI_Allreduce(&Parameters->neq, &nEq, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	/*************************************************************/
	//		Paraview output to file
	/************************************************************/
	Paraview_Output(RANK, NPROC, Parameters, FemStructs, FemFunctions);
	Paraview_Output_3D(RANK, NPROC, Parameters, FemStructs, FemFunctions);
	/*************************************************************/


	/****************************************************************************************/
	// 			Printing final result
	/****************************************************************************************/
	if (RANK == 0){
		sprintf(FileName,"../03_output/%02d/%s_%s_%s_%s_%s_%s_%s_N%d_E%d.txt", NPROC,Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->ShockCapture,
			Parameters->TimeIntegration, Parameters->MatrixVectorProductScheme, Parameters->Preconditioner, Parameters->nNodes, Parameters->nEl);
		OutFile = myfopen(FileName,"w");
		fprintf(OutFile, "\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
		fprintf(OutFile, "Problem Title: %s\n", Parameters->ProblemTitle);
		fprintf(OutFile, "Number of nodes: %d\n", Parameters->nNodes);
		fprintf(OutFile, "Number of elements: %d\n", Parameters->nEl);
		fprintf(OutFile, "Number of equations: %d\n", nEq);
		fprintf(OutFile, "Stabilization form used: %s\n", Parameters->StabilizationForm);
		fprintf(OutFile, "Discontinuities capture operator used: %s\n", Parameters->ShockCapture);
		fprintf(OutFile, "Stabilization coefficient tolerance: %E\n", Parameters->StabilizationTolerance);
		fprintf(OutFile, "Time integration method: %s\n", Parameters->TimeIntegration);
		fprintf(OutFile, "Time integration tolerance: %E\n", Parameters->TimeIntegrationTolerance);
		fprintf(OutFile, "Non linear tolerance: %E\n", Parameters->NonLinearTolerance);
		fprintf(OutFile, "Correction stopped: %s\n", Parameters->StopMulticorrection);
		fprintf(OutFile, "Maximum number of correction iteration: %d\n", Parameters->NonLinearMaxIter);
		fprintf(OutFile, "Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
		fprintf(OutFile, "Reordering: %s\n", Parameters->reordering);
		fprintf(OutFile, "Solver used: %s\n", Parameters->Solver);
		fprintf(OutFile, "Solver tolerance used: %E\n", Parameters->SolverTolerance);
		fprintf(OutFile, "Maximum number of solver iteration: %d\n", Parameters->LinearMaxIter);
		fprintf(OutFile, "Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
		fprintf(OutFile, "Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);
		fprintf(OutFile, "Number of MPI ranks: %d\n", NPROC);
		fprintf(OutFile, "Alpha: %lf\t Step time: %E\t Final Time: %lf - Stopped at Steady State? %s (Current Time: %lf)\n", Parameters->Alpha, Parameters->DeltaT,
			Parameters->FinalTime, Parameters->StopAtSteadyState, Parameters->CurrentTime);
		fprintf(OutFile, "\n========================================================================\n\n");
		fclose(OutFile);

		printf("\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
		printf("Problem Title: %s\n", Parameters->ProblemTitle);
		printf("Number of nodes: %d\n", Parameters->nNodes);
		printf("Number of elements: %d\n", Parameters->nEl);
		printf("Number of equations: %d\n", nEq);
		printf("Stabilization form used: %s\n", Parameters->StabilizationForm);
		printf("Discontinuities capture operator used: %s\n", Parameters->ShockCapture);
		printf("Stabilization coefficient tolerance: %E\n", Parameters->StabilizationTolerance);
		printf("Time integration method: %s\n", Parameters->TimeIntegration);
		printf("Time integration tolerance: %E\n", Parameters->TimeIntegrationTolerance);
		printf("Non linear tolerance: %E\n", Parameters->NonLinearTolerance);
		printf("Correction stopped: %s\n", Parameters->StopMulticorrection);
		printf("Maximum number of correction iteration: %d\n", Parameters->NonLinearMaxIter);
		printf("Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
		printf("Solver used: %s\n", Parameters->Solver);
		printf("Solver tolerance used: %E\n", Parameters->SolverTolerance);
		printf("Preconditioner used: %s\n", Parameters->Preconditioner);
		printf("Scaling used: %s\n", Parameters->Scaling);
		printf("Maximum number of solver iteration: %d\n", Parameters->LinearMaxIter);
		printf("Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
		printf("Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);
		printf("Number of MPI ranks: %d\n", NPROC);
		printf("Alpha: %lf\t Step time: %E\t Final Time: %lf - Stopped at Steady State? %s (Current Time: %lf)\n", Parameters->Alpha, Parameters->DeltaT,
		Parameters->FinalTime, Parameters->StopAtSteadyState, Parameters->CurrentTime);
		printf("\n========================================================================\n\n");
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else{
		MPI_Barrier(MPI_COMM_WORLD);
	}
	printf("Reordering: %s (bandwidth before: %d) (bandwidth after: %d) RANK: %d\n", Parameters->reordering, 
		Parameters->bandwidth_bef, Parameters->bandwidth_aft, RANK);

	/***************************************************************************************/
	//				Memory deallocation
	/**************************************************************************************/
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3) == 0){
		free(MatrixData->A);
		free(MatrixData->Aaux);
		free(MatrixData->Perm);
		free(MatrixData->invPerm);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
		if (strcasecmp(Parameters->Preconditioner,"BlockDiag")==0){
			free(MatrixData->Id);
			free(MatrixData->IdAux);
			free(MatrixData->BlockDiag);
			free(MatrixData->BlockDiagAux);
			free(MatrixData->invBlockDiag);
			free(MatrixData->invBlockDiagAux);
		}
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") == 0){
		int I, nel;

		nel = Parameters->nel;
		free(MatrixData->AA);
		free(MatrixData->JA);
		free(MatrixData->IA);
		free(MatrixData->Diag);
		free(MatrixData->invDiag);
		for (I = 0; I < nel; I++){
			free(MatrixData->Scheme_by_Element[I]);
			free(MatrixData->Scheme_by_Element_bef[I]);
			free(MatrixData->Scheme_by_Element_aft[I]);
		}
		free(MatrixData->Scheme_by_Element);
		free(MatrixData->Scheme_by_Element_bef);
		free(MatrixData->Scheme_by_Element_aft);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
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

	free(MatrixData);
	free(FemStructs->Node);
	free(FemStructs->Element);
	free(FemStructs->F);
	free(FemStructs->U);
	free(FemStructs->eqrho);
	if (RANK != 0){
		free(FemStructs->SendBuffer_bef[0]);
		free(FemStructs->SendBuffer_bef[1]);
		free(FemStructs->RecvBuffer_bef[0]);
		free(FemStructs->RecvBuffer_bef[1]);
		free(FemStructs->IdSend_bef);    	
		free(FemStructs->IdRecv_bef);
	}
	if (RANK != NPROC-1){
		free(FemStructs->SendBuffer_aft[0]);
		free(FemStructs->SendBuffer_aft[1]);
		free(FemStructs->RecvBuffer_aft[0]);
		free(FemStructs->RecvBuffer_aft[1]);
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


