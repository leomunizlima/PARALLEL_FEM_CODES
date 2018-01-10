#include "time_integration.h"

int Predictor_Old_TRBDF2(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
		FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int I, i;
	int nel, neq, passo;
	double t, dt, alpha, norm_a, norm_Da;
	double *a, *auxVec, *Da, *u, *u_old, *R; //Parametros do Preditor
	double *uB_old, *delta_old;
	double lambda, A, B, dtTR, dtBDF2;
	AuxBuildStructuresType *AuxBuild;

	nel = Parameters->nel;
	neq = Parameters->neq;
	
	a = (double*) mycalloc("a of 'Predictor_Old_TRBDF2'", neq + 1, sizeof(double));
	Da = (double*) mycalloc("Da of 'Predictor_Old_TRBDF2'", neq + 1, sizeof(double));
	auxVec = (double*) mycalloc("auxVec of 'Predictor_Old_TRBDF2'", neq + 1, sizeof(double));
	u_old = (double*) mycalloc("u_old of 'Predictor_Old_TRBDF2'", neq + 1, sizeof(double));
	R = (double*) mycalloc("R of 'Preditor_Old'", neq + 1, sizeof(double));
	delta_old = (double*) mycalloc("delta_old of 'Predictor_Old_TRBDF2'", nel, sizeof(double));
	uB_old = (double*) mycalloc("uB_old of 'Preditor_Old_TRBDF2'", nel*NDOF, sizeof(double));
	
	dt = Parameters->DeltaT;
	AuxBuild = (AuxBuildStructuresType*) mycalloc("AuxBuild of 'Predictor_Old_TRBDF2'",1,sizeof(AuxBuildStructuresType));
	AuxBuild->tolerance = Parameters->StabilizationTolerance;
	u = FemStructs->u;
	FemStructs->du = a;
	FemStructs->uB = uB_old;
	FemStructs->duB = NULL;
	FemStructs->delta_old = delta_old;

	//TRBDF2 parameters
	lambda = 2.0 - sqrt(2.0);
	A = 1.0/( lambda*( 2.0 - lambda ));
	B = ( 1.0 - lambda )*( 1.0 - lambda )/( lambda*( 2.0 - lambda ));

	dtTR = lambda*dt;
	dtBDF2 = (1.0 - lambda)*dt; // dt = dtTR + dtBDF2
	
	FemFunctions->InitialSolution(Parameters, FemStructs->Node, u);

	t = 0.0;
	passo = 0;
	int tag = 1;

	dmemcpy(neq, u, u_old); // copy u to u_old		
	Parameters->DeltaT_Build = dtTR;

	do{
		passo++;
		t += dtTR;
		#ifdef debug
			printf("\n\n Passo: %d\n", passo); 
		#endif
		/*-*-*-*-*-*-*-*- TR method -*-*-*-*-*-*-*-*/
		//PREDICAO
		Parameters->Alpha = 0.5;    //1.0/2.0;
		alpha = 0.5;     //1.0/2.0;
		Parameters->DeltaT = dtTR;
		i = 0;

		daxpy(neq, alpha*dtTR , a, u); // u  = u + alpha*dtTR*a
		memsetzero(neq,a); // set 0 in a
		
		norm_a = sqrt(ddot(neq, a, a));
		norm_Da = norm_a + 1.0;
		
		
		// MULTICORRECAO
		#ifdef debug
			printf("\nTR time integration:\n");				
		#endif
		do{
			i++;

			FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);

			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, R);

			FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, R, Da);

			daxpy(neq, 1, Da, a);
			daxpy(neq, alpha*dtTR, Da, u);
			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = sqrt(ddot(neq, Da, Da));
			#ifdef debug
				double normR;
				normR = sqrt(ddot(neq, R, R));
				printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n", 
					Parameters->NonLinearTolerance*norm_a, normR, norm_a, norm_Da, t, i);
			#endif
			}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection
			
		
			/*-*-*-*-*-*-*-*- BDF2 method -*-*-*-*-*-*-*-*/		
			
			//PREDICAO
			Parameters->Alpha = 1.0/( 2.0 - lambda );
			alpha = 1.0/( 2.0 - lambda );
			Parameters->DeltaT = dtBDF2;
			Parameters->DeltaT_Build = dtBDF2;
			t += dtBDF2;
			
			i = 0;

			for(I = 0; I < neq; I++)
			{
				//u_old[I] = u[I];
				u[I] = A*u[I] - B*u_old[I];// + alpha*dtBDF2*a[I]; 
				a[I] = 0.0;
			}
				
			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = norm_a + 1.0;
			
			// MULTICORRECAO
			#ifdef debug
				printf("\nBDF2 time integration:\n");				
			#endif
			do{
				i++;

				FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);

				FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, R);

				FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, R, Da);
				
				daxpy(neq, 1, Da, a);
				daxpy(neq, alpha*dtBDF2, Da, u);
				
				norm_a = sqrt(ddot(neq, a, a));
				norm_Da = sqrt(ddot(neq, Da, Da));
				#ifdef debug					
					double normR;
					normR = sqrt(ddot(neq, R, R));
					printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n", 
						Parameters->NonLinearTolerance*norm_a, normR, norm_a, norm_Da, t, i);
				#endif
			}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection
	
		dmemcpy(neq, u, u_old); // copy u to u_old		
	
		#ifdef debug
			printf("\n\n"); 
		#endif

	
	}while(!FemFunctions->StopTimeIntegration(Parameters,u,u_old,t)); // end while time

	free(a);
	free(Da);
	free(u_old);
	free(R);
	free(auxVec);
	free(uB_old);
	free(delta_old);
	free(AuxBuild);

	return 0;

}



