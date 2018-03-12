#include "NavierStokesEquations.h"
#include "../../../00_CommonFiles/CPU_Time_Operations/CPU_time.h"

int main(int argc, char **argv)
{
	ParametersType *Parameters=NULL;
	MatrixDataType *MatrixData;
	FemStructsType *FemStructs;
	FemFunctionsType *FemFunctions;
	FemOtherFunctionsType *FemOtherFunctions;
	struct timespec Start, End;
	double Preprocess_Time, Process_Time, Postprocess_Time;

	//*****************************//
	//	MPI Variables
	int RANK; // ID of MPI rank
	int NPROC; // Number of MPI ranks
	//*****************************//


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
	MPI_Comm_size(MPI_COMM_WORLD, &NPROC);

	/******************************************************** Preprocess *******************************************************/
	clock_gettime(CLOCK_REALTIME, &Start);
	Preprocess(RANK, NPROC, argc, argv, &Parameters, &MatrixData, &FemStructs, &FemFunctions, &FemOtherFunctions);
	clock_gettime(CLOCK_REALTIME, &End);
	Preprocess_Time = End.tv_sec - Start.tv_sec + 1e-9*(End.tv_nsec - Start.tv_nsec);

	/***************************************************************************************************************************/


	/********************************************************* Process *********************************************************/
	clock_gettime(CLOCK_REALTIME, &Start);
	Process(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
	clock_gettime(CLOCK_REALTIME, &End);
	Process_Time = End.tv_sec - Start.tv_sec + 1e-9*(End.tv_nsec - Start.tv_nsec);
	/***************************************************************************************************************************/

	
	/******************************************************** Postprocess *******************************************************/
	clock_gettime(CLOCK_REALTIME, &Start);
	Postprocess(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
	clock_gettime(CLOCK_REALTIME, &End);
	Postprocess_Time = End.tv_sec - Start.tv_sec + 1e-9*(End.tv_nsec - Start.tv_nsec);
	/***************************************************************************************************************************/

	
	/********************************************** Total Time ****************************************************************/
	calculateTime(Preprocess_Time, Process_Time, Postprocess_Time, Parameters);
	/**************************************************************************************************************************/

	free(Parameters);

	MPI_Finalize();

	return 0;
}


