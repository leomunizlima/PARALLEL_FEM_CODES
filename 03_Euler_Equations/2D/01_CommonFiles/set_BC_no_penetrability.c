#include "EulerEquations.h"
#include "../CYLINDER/cylinder.h"
#include "../NACA0012/naca0012.h"

void set_BC_no_penetrability(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->ProblemTitle,"CYLINDER")==0){
		FemFunctions->BC_theta = CYLINDER_theta;
		FemFunctions->BC_General_theta = BC_theta_OK;
		FemFunctions->BC_no_penetrability = CYLINDER_BC_no_penetrability;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"NACA0012")==0){
		FemFunctions->BC_theta = NACA0012_theta;
		FemFunctions->BC_General_theta = BC_theta_OK;
		FemFunctions->BC_no_penetrability = NACA0012_BC_no_penetrability;
	}
	else{
		FemFunctions->BC_theta = NULL;
		FemFunctions->BC_General_theta = BC_theta_NO;
		FemFunctions->BC_no_penetrability = NO_BC_no_penetrability;
	}
}


