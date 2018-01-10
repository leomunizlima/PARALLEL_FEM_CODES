#include "time_integration.h"

int Predictor_Old_BDF(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
		FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)

{
	int I, i;
	int nel, neq, passo;
	double t, dt, alpha, norm_a, norm_Da;
	double *a, *R, *Da,*u, *u_old, *u1, *u2; 
	double *uB_old, *delta_old, *delta_old_NMV;
	AuxBuildStructuresType *AuxBuild;

	nel = Parameters->nel;
	neq = Parameters->neq;
	
	u_old = (double*) mycalloc("u_old of 'Predictor_Old_BDF'", neq + 1, sizeof(double));
	a = (double*) mycalloc("a of 'Predictor_Old_BDF'", neq + 1, sizeof(double));
	Da = (double*) mycalloc("Da of 'Predictor_Old_BDF'", neq + 1, sizeof(double));	
	u1 = (double*) mycalloc("u1 of 'Predictor_Old_BDF'", neq + 1, sizeof(double));
	u2 = (double*) mycalloc("u2 of 'Predictor_Old_BDF'", neq + 1, sizeof(double));
	uB_old = (double*) mycalloc("uB_old of 'Predictor_Old_BDF'", nel*NDOF, sizeof(double));
	delta_old = (double*) mycalloc("delta_old of 'Predictor_BDF'", nel, sizeof(double));
	delta_old_NMV = (double*) mycalloc("delta_old of 'Preditor_Old_BDF'", nel, sizeof(double));

	u = FemStructs->u;
	R = FemStructs->F;
	dt = Parameters->DeltaT;
	Parameters->DeltaT_Build = dt;
	alpha = Parameters->Alpha;
	Parameters->Alpha_Build = alpha;
	AuxBuild = (AuxBuildStructuresType*) mycalloc("AuxBuild of 'Predictor_Old_BDF'",1,sizeof(AuxBuildStructuresType));
	AuxBuild->tolerance = Parameters->StabilizationTolerance;
	AuxBuild->delta_old_NMV = delta_old_NMV;
	FemStructs->AuxBuild = AuxBuild;
	FemStructs->du = a;
	FemStructs->uB = uB_old;
	FemStructs->duB = NULL;
	FemStructs->delta_old = delta_old;

	FemFunctions->InitialSolution(Parameters, FemStructs->Node, u);

	t = 0.0;
	passo = 0;
	int tag = 1;


	/*-----------------------------------------------------------------*/
	//                              BDF1
	/*-----------------------------------------------------------------*/
	
	passo++;
	t += dt;
	Parameters->Alpha = 1.0;	

	#ifdef debug
		printf("\n\n Passo: %d\n", passo); 
	#endif
	
	//PREDICAO
	i = 0;
	

	dmemcpy(neq, u, u_old); // copy u to u_old		
	dmemcpy(neq, u, u1); // copy u to u1		
	memsetzero(neq, a); // set 0 in a

 
	// MULTICORRECAO
	do{
		i++;
		
		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
		
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, R);

		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, R, Da);

		daxpy(neq, 1, Da, a);			//a^{i+1} = a^{i} + Da		-no mesmo passo de tempo
		daxpy(neq, alpha*dt, Da, u);            //u^{i+1} = u^{i} + dt*Da 	-no mesmo passo de tempo
		
		norm_a = sqrt(ddot(neq, a, a));
		norm_Da = sqrt(ddot(neq, Da, Da));

		#ifdef debug
			double normR;
			normR = sqrt(ddot(neq, R, R));
			printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n",
					Parameters->NonLinearTolerance*norm_a, normR, norm_a, norm_Da, t, i);
		#endif
	}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection
	
	dmemcpy(neq, u, u2); // copy u to u2		
	
	/*-----------------------------------------------------------------*/
	//                              BDF2
	/*-----------------------------------------------------------------*/
	double four_thirds = 4.0/3.0;
	double two_thirds = 2.0/3.0;
	double one_third = 1.0/3.0;
	Parameters->Alpha = two_thirds; 

	do{
		passo++;
		t += dt;

		#ifdef debug
			printf("\n\n Passo: %d\n", passo); 
		#endif
		
		//PREDICAO
		i = 0;
		
		for(I = 0; I < neq; I++)
		{
			u[I] = four_thirds*u2[I] - one_third*u1[I];
			a[I] = 0.0;
		}
		 
		// MULTICORRECAO
		do{
			i++;
			
			FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
			
			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, R);

			FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, R, Da);
			
			daxpy(neq, 1, Da, a);          		//a^{i+1} = a^{i} + Da			-no mesmo passo de tempo
			daxpy(neq, two_thirds*dt, Da, u);          //u^{i+1} = u^{i} + 2/3dt*Da 	-no mesmo passo de tempo
			
			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = sqrt(ddot(neq, Da, Da));
			#ifdef debug
				double normR;
				normR = sqrt(ddot(neq, R, R));
				printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n",
						Parameters->NonLinearTolerance*norm_a, normR, norm_a, norm_Da, t, i);
			#endif
		}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection
	
		dmemcpy(neq, u2, u1); // copy u2 to u1		
		dmemcpy(neq, u, u2); // copy u to u2		

		#ifdef debug
			printf("\n\n"); 
		#endif

	}while(!FemFunctions->StopTimeIntegration(Parameters,u,u_old,t)); // end while time

	free(u_old);
	free(a);
	free(Da);
	free(u1);
	free(u2);	
	free(uB_old);
	free(delta_old);
	free(delta_old_NMV);
	free(AuxBuild);

	return 0;

}



