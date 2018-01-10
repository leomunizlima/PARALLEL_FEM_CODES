#include "time_integration.h"

int Predictor_New_BDF(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
		FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int I, i;
	int nel, neq, passo;
	double t, dt, alpha, norm_a, norm_Da;
	double *a, *aB, *Da, *DaB, *u, *u_old, *u1, *u2, **R2, *R2aux, *invN2, **M2, *M2aux, *F; //Parametros do Preditor
	double *uB, *uB1, *uB2, *delta_old, *delta_old_NMV;

	AuxBuildStructuresType *AuxBuild;
	
	nel = Parameters->nel;
	neq = Parameters->neq;
	
	u_old = (double*) mycalloc("u_old of 'Preditor_New'", neq + 1, sizeof(double));
	a = (double*) mycalloc("a of 'Preditor_New'", neq + 1, sizeof(double));
	aB = (double*) mycalloc("aB of 'Preditor_New'",NDOF*nel, sizeof(double));
	Da = (double*) mycalloc("Da of 'Preditor_New'", neq + 1, sizeof(double));
	DaB = (double*) mycalloc("DaB of 'Preditor_New'", NDOF*nel, sizeof(double));
	u1 = (double*) mycalloc("u1 of 'Preditor_New'", neq + 1, sizeof(double));
	u2 = (double*) mycalloc("u2 of 'Preditor_New'", neq + 1, sizeof(double));
	R2 = (double**) mycalloc("R2 of 'Preditor_New'", nel, sizeof(double));
	R2aux = (double *) mycalloc("R2aux of 'Preditor_New'", NDOF*nel, sizeof(double));
	for (I = 0; I < nel;I++)
		R2[I] = &R2aux[NDOF*I]; 
	invN2 = (double*) mycalloc("invN2 of 'Preditor_New'", nel, sizeof(double));
	M2 = (double **) mycalloc("M2 of 'Preditor_New'", nel, sizeof(double*));
	M2aux = (double *) mycalloc("M2aux of 'Preditor_New'", NNOEL*NDOF*NDOF*nel, sizeof(double));
	for (I = 0; I < nel;I++)
		M2[I] = &M2aux[NNOEL*NDOF*NDOF*I]; 
	uB = (double*) mycalloc("uB of 'Preditor_New'", nel*NDOF, sizeof(double));
	uB1 = (double*) mycalloc("uB1 of 'Preditor_New'", nel*NDOF, sizeof(double));
	uB2 = (double*) mycalloc("uB2 of 'Preditor_New'", nel*NDOF, sizeof(double));
	delta_old = (double*) mycalloc("delta_old of 'Preditor_New'", nel, sizeof(double));
	delta_old_NMV = (double*) mycalloc("delta_old of 'Preditor_New'", nel, sizeof(double));

	u = FemStructs->u;
	F = FemStructs->F;
	dt = Parameters->DeltaT;
	Parameters->DeltaT_Build = dt;
	AuxBuild = (AuxBuildStructuresType*) mycalloc("AuxBuild of 'Predictor_New_BDF'",1,sizeof(AuxBuildStructuresType));
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
	
	FemFunctions->InitialSolution(Parameters, FemStructs->Node, u);
	//uB_InitialSolution(Parameters, FemStructs, FemFunctions, u, uB);
	
	t = 0.0;
	passo = 1;
	int tag = 1;

	t += dt;

	#ifdef debug
		printf("\n\n Passo: %d\n", passo); 
	#endif

	//PREDICAO
	Parameters->Alpha = 1.0;
	alpha = 1.0;
	Parameters->Alpha_Build = alpha;
	i = 0;

	dmemcpy(neq, u, u1); // copy u to u1		
	memsetzero(neq, a); // set 0 in a

	dmemcpy(nel*NDOF, uB, uB1); // copy uB to uB1		
	memsetzero(nel*NDOF, aB); // set 0 in a
		
	// MULTICORRECAO
	do{
		i++;

		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
		
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, F);

		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, Da);

		calculate_DaB(Parameters, FemStructs, Da, DaB);
		
		daxpy(neq, 1, Da, a);           //a^{i+1} = a^{i} + Da		-no mesmo passo de tempo
		daxpy(neq, dt, Da, u);          //u^{i+1} = u^{i} + dt*Da 	-no mesmo passo de tempo
		
		daxpy(nel*NDOF, 1, DaB, aB);           //aB^{i+1} = aB^{i} + DaB		-no mesmo passo de tempo
		daxpy(nel*NDOF, dt, DaB, uB);          //uB^{i+1} = uB^{i} + dt*DaB 	-no mesmo passo de tempo
	
		norm_a = sqrt(ddot(neq, a, a));
		norm_Da = sqrt(ddot(neq, Da, Da));

		#ifdef debug
			double normF;
			normF = sqrt(ddot(neq, F, F));
			printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n", 
				Parameters->NonLinearTolerance*norm_a, normF, norm_a, norm_Da, t, i);
		#endif
	}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection

	dmemcpy(neq, u, u2); // copy u to u2		
	dmemcpy(nel*NDOF, uB, uB2); // copy uB to uB2		

	/*-----------------------------------------------------------------*/
	//                              BDF2
	/*-----------------------------------------------------------------*/
	double four_thirds = 4.0/3.0;
	double two_thirds = 2.0/3.0;
	double one_third = 1.0/3.0;
	do{
		passo++;
		t += dt;
	
		#ifdef debug
			printf("\n\n Passo: %d\n", passo); 
		#endif

		//PREDICAO
		Parameters->Alpha = 2.0/3.0;
		alpha = 2.0/3.0;
		Parameters->Alpha_Build = alpha;
		i = 0;

		dmemcpy(neq, u, u_old); // copy u to u_old		

		for(I = 0; I < neq; I++)
		{
			u[I] = four_thirds*u2[I] - one_third*u1[I]; 
			a[I] = 0.0;
		}
			
		for(I = 0; I < nel*NDOF; I++){
			uB[I] = four_thirds*uB2[I] - one_third*uB1[I];
			aB[I] = 0.0;
		}
		 
			
		// MULTICORRECAO
		do{
			i++;

			FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);

			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, F);

			FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, Da);

			calculate_DaB(Parameters, FemStructs, Da, DaB);

			daxpy(neq, 1, Da, a);
			daxpy(neq, two_thirds*dt, Da, u);

			daxpy(nel*NDOF, 1, DaB, aB);
			daxpy(nel*NDOF, two_thirds*dt, DaB, uB);

			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = sqrt(ddot(neq, Da, Da));
			#ifdef debug
				double normF;
				normF = sqrt(ddot(neq, F, F));
				printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n", 
					Parameters->NonLinearTolerance*norm_a, normF, norm_a, norm_Da, t, i);
			#endif
		}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection

		dmemcpy(neq, u2, u1); // copy u2 to u1		
		dmemcpy(neq, u, u2); // copy u to u2		
	
		dmemcpy(nel*NDOF, uB2, uB1); // copy uB2 to uB1		
		dmemcpy(nel*NDOF, uB, uB2); // copy uB to uB2		
	
		#ifdef debug
			printf("\n\n"); 
		#endif

	}while(!FemFunctions->StopTimeIntegration(Parameters,u,u_old,t)); // end while time

	free(u_old);
	free(a);
	free(Da);
	free(uB);
	free(uB1);
	free(uB2);
	free(u1);
	free(u2);
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



