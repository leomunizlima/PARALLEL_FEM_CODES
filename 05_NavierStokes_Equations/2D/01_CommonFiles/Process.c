#include "NavierStokesEquations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int i;
	int neq_bef = Parameters->neq_bef; 
	int neq = Parameters->neq;
	int RANK, NPROC;
	int NonLinearMaxIter = Parameters->NonLinearMaxIter;
	double NonLinearTolerance = Parameters->NonLinearTolerance; 
	double norms, normu, delta2;	
	double *u = FemStructs->u;
	double *F = FemStructs->F;
	double *s;

	setProblem(Parameters, FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters,FemOtherFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions);
	RANK = Parameters->RANK;
	NPROC = Parameters->NPROC;

	//=========Iteração de Newton===========
	s = (double*) mycalloc("s of 'Process'", Parameters->NEQ+1, sizeof(double));

	i = 0;
	delta2 = NonLinearTolerance + 1.0;

	while(delta2 > NonLinearTolerance && i < NonLinearMaxIter){
		i++;		
		MPI_VectorUpdate(RANK, NPROC, Parameters, FemStructs, u, 0);
		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, i, F);
		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, s);
		daxpy(neq, 1, &s[neq_bef], u);		//u = uold + s
		
		//====== Condicoes de saida ======		

		norms = sqrt(ddot(neq, &s[neq_bef], &s[neq_bef]));
		normu = sqrt(ddot(neq, u, u));
		double normF = sqrt(ddot(neq,&F[neq_bef],&F[neq_bef]));
		delta2 = norms / normu;
		if (RANK==0){ 
			printf("===================================================");	
			printf("\n         |s|/|u| = %3.2E |s|=%lf |u|=%lf |F|=%lf (iteration %d)\n\n", delta2, norms, normu, normF, i);
		}
					
	}

	if (RANK==0){ 
		printf("===================================================\n");
		printf("\n Newton iterations = %d  \n", i);
	}
		
	free(s);
	
	return 0;
}


















