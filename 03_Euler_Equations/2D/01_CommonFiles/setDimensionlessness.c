#include "EulerEquations.h"

void setDimensionlessness(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Dimensionless,"YES")==0)
		FemFunctions->Ax_Ay_calculations = dimensionless_Ax_Ay_calculations;
	else if (strcasecmp(Parameters->Dimensionless,"NO")==0)
		FemFunctions->Ax_Ay_calculations = dimensional_Ax_Ay_calculations;
	else{
		printf("Dimensionless is not defined correctly!\n");
		exit(1);
	}	

}
