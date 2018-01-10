#include "CPU_time.h"

void calculateTime(double Preprocess_Time, double Process_Time, double Postprocess_Time, ParametersType *Parameters)
{
	int h, m;
	double s, Total_Time;
	FILE *OutFile;
	char FileName[300];
	int RANK, NPROC;

	RANK = Parameters->RANK;
	NPROC = Parameters->NPROC;

	MPI_Barrier(MPI_COMM_WORLD);
	if (RANK == 0){
		Total_Time = Preprocess_Time + Process_Time + Postprocess_Time;
		h = (int) Total_Time;
		h = h/3600;
		m = (int) (Total_Time - 3600*h);
		m = m/60;
		s = Total_Time - (3600*h + 60*m);
        
		#ifdef SSTranspEquation2D
			sprintf(FileName, "../03_output/NPROC_%d_%s_%s_%s_%s_%s_%s_%s_N%d_E%d.dat",NPROC, Parameters->ProblemTitle,Parameters->StabilizationForm,
				Parameters->ShockCapture, Parameters->h_Shock, Parameters->MatrixVectorProductScheme,Parameters->Solver,Parameters->Preconditioner, 
				Parameters->nNodes,Parameters->nEl); 	
		#endif
        
		#ifdef TranspEquation2D
			sprintf(FileName, "../03_output/NPROC_%d_%s_%s_%s_%s_%s_%s_%s_%s_N%d_E%d.dat",NPROC, Parameters->ProblemTitle,Parameters->TimeIntegration,
				Parameters->StabilizationForm, Parameters->ShockCapture, Parameters->h_Shock, Parameters->MatrixVectorProductScheme,Parameters->Solver,
				Parameters->Preconditioner,Parameters->nNodes,Parameters->nEl); 	
		#endif
        
		#ifdef EulerEquations2D
			sprintf(FileName,"../03_output/NPROC_%d_%s_%s_%s_%s_%s_%s_%s_N%d_E%d.dat", NPROC, Parameters->Experiments, Parameters->ProblemTitle, 
				Parameters->StabilizationForm, Parameters->ShockCapture, Parameters->TimeIntegration, Parameters->MatrixVectorProductScheme,
				Parameters->Preconditioner,Parameters->nNodes, Parameters->nEl);
		#endif
		
		#ifdef SSNavierStokesEquations2D
			sprintf(FileName,"../03_output/NPROC_%d_%s_%s_%s_%s_%s_N%d_E%d.dat", NPROC, Parameters->Experiments, Parameters->ProblemTitle, 
				Parameters->StabilizationForm, Parameters->MatrixVectorProductScheme,Parameters->Preconditioner,Parameters->nNodes, Parameters->nEl);
		#endif

		OutFile = myfopen(FileName,"a");
		fprintf(OutFile,"\n============================ Processing Time ============================\n\n");
		fprintf(OutFile,"\nPreprocess Time: %lf\n", Preprocess_Time);	
		fprintf(OutFile,"\nProcess Time: %lf\n", Process_Time);	
		fprintf(OutFile,"\nPostprocess Time: %lf\n", Postprocess_Time);	
		fprintf(OutFile,"\nTOTAL TIME: %lf (%dh %dm %lfs)\n", Total_Time, h, m, s);
		fprintf(OutFile,"\n=========================================================================\n\n");
		fclose(OutFile);
		
		printf("\n============================ Processing Time ============================\n\n");
		printf("\nPreprocess Time: %lf (RANK: %d)\n", Preprocess_Time, RANK);	
		printf("\nProcess Time: %lf (RANK: %d)\n", Process_Time, RANK);	
		printf("\nPostprocess Time: %lf (RANK: %d)\n", Postprocess_Time, RANK);	
		printf("\nTOTAL TIME: %lf (%dh %dm %lfs) (RANK: %d)\n", Total_Time, h, m, s, RANK);
		printf("\n=========================================================================\n\n");
	}	
}



