#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/preconditioners.h"

int Process(int RANK, int NPROC, ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int i, neq = Parameters->neq;
	int nel = Parameters->nel;
	double NonLinearTolerance = Parameters->NonLinearTolerance; 
	double normdiff, normu;	
	double *u = FemStructs->u;
	double *U = FemStructs->U;
	double *F = FemStructs->F;
	double *uold = mycalloc("uold",sizeof(double),neq+1);	
	double *delta_old = mycalloc("delta_old",sizeof(double),nel);
	double *Re_old = mycalloc("Re_old",sizeof(double),nel);

	setProblem(Parameters, FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters, FemOtherFunctions);
	setScaling(Parameters, FemFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions);
	
	FemStructs->delta_old = delta_old;
	FemStructs->Re_old = Re_old;

	normdiff = NonLinearTolerance + 1;
	normu = 1.0/NonLinearTolerance;
	i = 0;

/*	int j;
	FILE *Out;
	char Filename[200];

	sprintf(Filename,"RANK_%d.dat",RANK);
	Out = fopen(Filename,"w");
*/
	while (normdiff > NonLinearTolerance*normu){
		i++;	
		dcopy(neq, u, uold);
		MPI_VectorUpdate(RANK, NPROC, Parameters, FemStructs, u, 0);
		FemOtherFunctions->Build_K_F(Parameters, MatrixData, FemStructs, FemFunctions);
		
/*		MPI_VectorUpdate(RANK, NPROC, Parameters, FemStructs, FemStructs->f, 0);

		fprintf(Out,"before:\n");
		for (i=0; i<Parameters->NEQ; i++)
			fprintf(Out,"F[%d]=%lf\n", i, FemStructs->F[i]);
		
		for (i=0; i<Parameters->nel;i++)
			fprintf(Out,"e:%d LM -> %d %d %d\n", i, FemStructs->lm[i][0], FemStructs->lm[i][1],FemStructs->lm[i][2]);	
		for (i=0; i<Parameters->nel;i++)
			for (j=0;j<9;j++)
				fprintf(Out,"A[%d][%d]=%lf\n",i,j,MatrixData->A[i][j]);				i=1;
*/
	
		FemFunctions->scaling(Parameters, MatrixData, FemStructs);
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, i, F);
		
		

/*		fprintf(Out,"after:\n");
		for (i=0; i<Parameters->NEQ; i++)
			fprintf(Out,"F[%d]=%lf\n", i, FemStructs->F[i]);		
		for (i=0; i<Parameters->nel;i++)
			for (j=0;j<9;j++)
				fprintf(Out,"LUe[%d]=%lf  A[%d][%d]=%lf\n",j,MatrixData->LUe[i][j],i,j,MatrixData->A[i][j]);			

		fclose(Out);
		
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,1);
*/
		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, U);
		FemFunctions->unscaling(Parameters, MatrixData, FemStructs, FemStructs->u);
		daxpy(neq, -1, u, uold);
		normdiff = sqrt(ddot(neq,uold,uold));
		normu = sqrt(ddot(neq,u,u));		

		#ifdef debug
			if (RANK==0) printf( "normdiff=%lf i=%d (normdiff/normu = %lf)\n\n\n", normdiff, i, normdiff/normu);
		#endif
	
	}

	free(Re_old);
	free(delta_old);	
	free(uold);

	return 0;
}


