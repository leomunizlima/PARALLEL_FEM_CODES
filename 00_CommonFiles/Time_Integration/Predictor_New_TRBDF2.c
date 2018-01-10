#include "time_integration.h"

int Predictor_New_TRBDF2(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
		FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int I, i;
	int nel, neq, passo;
	double t, dt, alpha, norm_a, norm_Da;
	double *a, *aB, *auxVec, *Da, *DaB, *u, *u_old, **R2, *R2aux, *invN2, **M2, *M2aux, *F; //Parametros do Preditor
	double *uB, *uB_old, *delta_old, *delta_old_NMV;
	double lambda, A, B, dtTR, dtBDF2;
	AuxBuildStructuresType *AuxBuild;

	nel = Parameters->nel;
	neq = Parameters->neq;
	
	a = (double*) mycalloc("a of 'Predictor_TRBDF2'", neq + 1, sizeof(double));
	aB = (double*) mycalloc("aB of 'Predictor_TRBDF2'",NDOF*nel, sizeof(double));
	Da = (double*) mycalloc("Da of 'Predictor_TRBDF2'", neq + 1, sizeof(double));
	DaB = (double*) mycalloc("DaB of 'Predictor_TRBDF2'", NDOF*nel, sizeof(double));
	auxVec = (double*) mycalloc("auxVec of 'Predictor_TRBDF2'", neq + 1, sizeof(double));
	u_old = (double*) mycalloc("u_old of 'Predictor_TRBDF2'", neq + 1, sizeof(double));
	R2 = (double**) mycalloc("R2 of 'Predictor_TRBDF2'", nel, sizeof(double));
	R2aux = (double *) mycalloc("R2aux of 'Predictor_TRBDF2'", NDOF*nel, sizeof(double));
	for (I = 0; I < nel;I++)
		R2[I] = &R2aux[NDOF*I]; 
	invN2 = (double*) mycalloc("invN2 of 'Predictor_TRBDF2'", nel, sizeof(double));
	M2 = (double **) mycalloc("M2 of 'Predictor_TRBDF2'", nel, sizeof(double*));
	M2aux = (double *) mycalloc("M2aux of 'Predictor_TRBDF2'", NNOEL*NDOF*NDOF*nel, sizeof(double));
	for (I = 0; I < nel;I++)
		M2[I] = &M2aux[NNOEL*NDOF*NDOF*I]; 
	uB = (double*) mycalloc("uB of 'Predictor_TRBDF2'", nel*NDOF, sizeof(double));
	uB_old = (double*) mycalloc("uB of 'Predictor_TRBDF2'", nel*NDOF, sizeof(double));
	delta_old = (double*) mycalloc("delta_old of 'Predictor_TRBDF2'", nel, sizeof(double));
	delta_old_NMV = (double*) mycalloc("delta_old of 'Preditor_New'", nel, sizeof(double));
	
	u = FemStructs->u;
	F = FemStructs->F;
	dt = Parameters->DeltaT;

	AuxBuild = (AuxBuildStructuresType*) mycalloc("AuxBuild of 'Predictor_New'",1,sizeof(AuxBuildStructuresType));
	AuxBuild->tolerance = Parameters->StabilizationTolerance;
	AuxBuild->M2 = M2;
	AuxBuild->R2 = R2;
	AuxBuild->invN2 = invN2;
	AuxBuild->delta_old_NMV = delta_old_NMV;
	FemStructs->AuxBuild = AuxBuild;
	FemStructs->delta_old = delta_old;
	FemStructs->du = a;
	FemStructs->uB = uB;
	FemStructs->duB = aB;
	

	//TRBDF2 parameters
	lambda = 2.0 - sqrt(2.0);
	A = 1.0/( lambda*( 2.0 - lambda ));
	B = ( 1.0 - lambda )*( 1.0 - lambda )/( lambda*( 2.0 - lambda ));

	dtTR = lambda*dt;
	dtBDF2 = (1.0 - lambda)*dt; // dt = dtTR + dtBDF2
	
	FemFunctions->InitialSolution(Parameters, FemStructs->Node, u);
	uB_InitialSolution(Parameters, FemStructs, FemFunctions, u, uB);

	t = 0.0;
	passo = 0;
	int tag = 1;

	dmemcpy(neq, u, u_old); // copy u to u_old		


	do{
		Parameters->DeltaT = dtTR;
		Parameters->DeltaT_Build = dtTR;
		Parameters->Alpha = 0.5;    //1.0/2.0;
		alpha = 0.5;     //1.0/2.0;
		Parameters->Alpha_Build = alpha;
		passo++;
		t += dtTR;
		#ifdef debug
			printf("\n\n Passo: %d\n", passo); 
		#endif
		/*-*-*-*-*-*-*-*- TR method -*-*-*-*-*-*-*-*/
		//PREDICAO

		i = 0;

		daxpy(neq, alpha*dtTR , a, u); // u  = u + alpha*dtTR*a
		memsetzero(neq,a); // set 0 in a
		
		daxpy(nel*NDOF, alpha*dtTR , aB, uB); // uB  = uB + alpha*dtTR*aB
		memsetzero(nel*NDOF,aB); // set 0 in aB
	
		norm_a = sqrt(ddot(neq, a, a));
		norm_Da = norm_a + 1.0;
		
		
		// MULTICORRECAO
		#ifdef debug
			printf("\nTR time integration:\n");				
		#endif


		do{
			i++;

			FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);

			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, F);

			FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, Da);

			calculate_DaB(Parameters, FemStructs, Da, DaB);

			daxpy(neq, 1, Da, a);
			daxpy(neq, alpha*dtTR, Da, u);

			daxpy(nel*NDOF, 1, DaB, aB);
			daxpy(nel*NDOF, alpha*dtTR, DaB, uB);

			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = sqrt(ddot(neq, Da, Da));
			#ifdef debug
				double normF;
				normF = sqrt(ddot(neq, F, F));
				printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n",
					Parameters->NonLinearTolerance*norm_a, normF, norm_a, norm_Da, t, i);
			#endif
		}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection
		
	
		/*-*-*-*-*-*-*-*- BDF2 method -*-*-*-*-*-*-*-*/		
		
		//PREDICAO
		Parameters->Alpha = 1.0/( 2.0 - lambda );
		alpha = 1.0/( 2.0 - lambda );
		Parameters->Alpha_Build = alpha;
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
			
		for(I = 0; I < nel*NDOF; I++){
			uB[I] = A*uB[I] - B*uB_old[I];// + alpha*dtBDF2*aB[I]; 
			aB[I] = 0.0;
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

			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, F);

			FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, Da);

			calculate_DaB(Parameters, FemStructs, Da, DaB);
			
			daxpy(neq, 1, Da, a);
			daxpy(neq, alpha*dtBDF2, Da, u);
			
			daxpy(nel*NDOF , 1, DaB, aB);
			daxpy(nel*NDOF, alpha*dtBDF2, DaB, uB);
			
			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = sqrt(ddot(neq, Da, Da));
			#ifdef debug					
				double normF;
				normF = sqrt(ddot(neq, F, F));
				printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n", 
					Parameters->NonLinearTolerance*norm_a, normF, norm_a, norm_Da, t, i);
			#endif
		}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection

		dmemcpy(neq, u, u_old); // copy u to u_old		
		dmemcpy(nel*NDOF, uB, uB_old); // copy uB to uB_old		

		#ifdef debug
			printf("\n\n"); 
		#endif
	
	}while(!FemFunctions->StopTimeIntegration(Parameters,u,u_old,t)); // end while time

	free(a);
	free(Da);
	free(u_old);
	free(auxVec);
	free(uB);
	free(uB_old);
	free(DaB);
	free(delta_old);
	free(M2aux);
	free(M2);
	free(R2aux);
	free(R2);
	free(invN2);
	free(delta_old_NMV);
	free(AuxBuild);
	
	return 0;

}





