#include "time_integration.h"

void setStopCriteria(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if(strcasecmp(Parameters->StopMulticorrection,"ITERATION") == 0)
		FemFunctions->StopCriteria = StopByIterations;
	else if(strcasecmp(Parameters->StopMulticorrection,"NORM") == 0)
		FemFunctions->StopCriteria = StopByNorm;
	else{
		printf("Stop criteria is not defined!\n");
		exit(1);
	}

	if(strcasecmp(Parameters->StopAtSteadyState,"YES") == 0)
		FemFunctions->StopTimeIntegration = StopBySteadyState;
	else if(strcasecmp(Parameters->StopAtSteadyState,"NO") == 0)
		FemFunctions->StopTimeIntegration = StopByTime;
	else{
		printf("Stop time integration is not defined!\n");
		exit(1);
	}

}

int StopByIterations(ParametersType *Parameters, double norm_a, double norm_Da, int i)
{
	if (i>=Parameters->NonLinearMaxIter)
		return 1;
	else
		return 0;
}



int StopByNorm(ParametersType *Parameters, double norm_a, double norm_Da, int i)
{
	if (norm_Da < (Parameters->NonLinearTolerance)*norm_a || i>=Parameters->NonLinearMaxIter)
		return 1;
	else
		return 0;
}


int StopBySteadyState(ParametersType *Parameters, double *u, double *u_old, double t)
{
	int i, neq;
	double norm_diff, norm_u;
	double *diff = mycalloc("diff of 'StopBySteadyState'",Parameters->neq,sizeof(double));
	
	neq = Parameters->neq;

	for (i=0; i<neq; i++)
		diff[i] = u[i] - u_old[i];

	norm_diff = sqrt(ddot(neq, diff, diff)); 
	norm_u = sqrt(ddot(neq,u,u));	

	free(diff);

	if ((Parameters->TimeIntegrationTolerance)*norm_u > norm_diff || fabs(Parameters->FinalTime-t) < Parameters->DeltaT){
		Parameters->CurrentTime = t;	
		return 1;
	}
	else
		return 0;			
			
}

int StopByTime(ParametersType *Parameters, double *u, double *u_old, double t)
{
	//REMOVER DEPOIS...---------------
	Parameters->CurrentTime = t;	
	//--------------------------------
	if (Parameters->FinalTime-t < 1e-10){
		Parameters->CurrentTime = t;	
		return 1;
	}
	else
		return 0;
}


