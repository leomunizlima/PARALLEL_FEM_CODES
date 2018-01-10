#include "time_integration.h"
////////////////////////////////////////////////
#include "../MatrixVector_Operations/matvec.h"
///////////////////////////////////////////////

int Predictor_Old(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
		FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)

{
	int i;
	int nel, neq, neq_bef, NEQ, passo;
	double t, dt, alpha, norm_a, norm_Da;
	double *a, *Da, *DAa, *u, *DU, *u_old, *R; //Parametros do Preditor
	double *uB_old, *delta_old;
	AuxBuildStructuresType *AuxBuild;

	int RANK = Parameters->RANK;
	int NPROC = Parameters->NPROC;

	nel = Parameters->nel;
	neq = Parameters->neq;
	neq_bef = Parameters->neq_bef;
	NEQ = Parameters->NEQ;

	DU = (double*) mycalloc("a of 'Preditor_Old'", NEQ + 1, sizeof(double));
	DAa = (double*) mycalloc("Da of 'Preditor_Old'", NEQ + 1, sizeof(double));
	u_old = (double*) mycalloc("u_old of 'Preditor_Old'", neq + 1, sizeof(double));
	uB_old = (double*) mycalloc("uB_old of 'Preditor_Old'", nel*NDOF, sizeof(double));
	delta_old = (double*) mycalloc("delta_old of 'Preditor_Old'", nel, sizeof(double));

	dt = Parameters->DeltaT;
	alpha = Parameters->Alpha;
	Parameters->DeltaT_Build = dt;
	Parameters->Alpha_Build = alpha;
	AuxBuild = (AuxBuildStructuresType*) mycalloc("AuxBuild of 'Predictor_Old'",1,sizeof(AuxBuildStructuresType));
	AuxBuild->tolerance = Parameters->StabilizationTolerance;
	FemStructs->AuxBuild = AuxBuild;
	u = FemStructs->u;
	R = FemStructs->F;
	Da = &DAa[neq_bef];
	a = &DU[neq_bef];
	FemStructs->du = a;
	FemStructs->uB = uB_old;
	FemStructs->duB = NULL;
	FemStructs->delta_old = delta_old;

	FemFunctions->InitialSolution(Parameters, FemStructs->Node, u);

	t = 0.0;
	passo = 0;
	int tag = 1;

	do{
		passo++;
		t += dt;

		#ifdef debug
			if (RANK == 0)
				printf("\n\n Passo: %d\n", passo); 
		#endif

		//PREDICAO
		i = 0;
		dcopy(neq, u, u_old);
		daxpy(neq, (1.0-alpha)*dt, a, u);
		dzero(neq, a);		
		// MULTICORRECAO
		do{
			i++;

			MPI_VectorUpdate(RANK, NPROC, Parameters, FemStructs, u, 0);
			MPI_VectorUpdate(RANK, NPROC, Parameters, FemStructs, a, 1);

			FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		
//		ebemvNDOF4(Parameters, MatrixData, FemStructs, FemStructs->U, DU);
//		ebemv(Parameters, MatrixData, FemStructs, FemStructs->U, DU);
		
	/*	double norm_u = sqrt(ddot(neq, u, u));
		double norm_R = sqrt(ddot(neq,&R[neq_bef],&R[neq_bef]));
		norm_a = sqrt(ddot(neq, a, a));
		if (RANK==0){
			printf("\n\n------\n");
			printf("norm_u=%.15lf\n",norm_u);
			printf("norm_R=%.15lf\n",norm_R);
			printf("norm_a=%.15lf\n",norm_a);
			printf("\n------\n\n");
		}*/
//		MPI_Abort(MPI_COMM_WORLD,1);	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, R);
			FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, R, DAa);


			daxpy(neq, 1, Da, a);
			daxpy(neq, alpha*dt, Da, u);
			
		
			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = sqrt(ddot(neq, Da, Da));
		
		//	MPI_VectorUpdate(RANK, NPROC, Parameters, FemStructs, u);
		
			#ifdef debug
				double normR;
				normR = sqrt(ddot(neq, &R[neq_bef], &R[neq_bef]));
				if (RANK == 0){
					printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n", 
					Parameters->NonLinearTolerance*norm_a, normR, norm_a, norm_Da, t, i);
				}
			#endif
		}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection
	
		#ifdef debug
			if (RANK == 0)
				printf("\n\n");
		#endif
	
	}while(!FemFunctions->StopTimeIntegration(Parameters,u,u_old,t)); // end while time

	free(DU);
	free(DAa);
	free(u_old);
	free(uB_old);
	free(delta_old);
	free(AuxBuild);

	return 0;

}



