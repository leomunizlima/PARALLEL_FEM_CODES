#include "EulerEquations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/scaling.h"

int setScaling(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Scaling,"NOT")==0){
		FemFunctions->unscaling = NO_unscaling;
		FemFunctions->scaling = NO_scaling;
	}
	else if (strcasecmp(Parameters->Scaling,"BlockScaling")==0){
		FemFunctions->unscaling = NO_unscaling;
		if (strcasecmp(Parameters->MatrixVectorProductScheme, "EBE")==0) 
			FemFunctions->scaling = Block_scaling_EBE;
		else{
			printf("Scaling is used only with EBE scheme!\n");
			exit(1);
		}
			
	}
	else{
		printf("Scaling is not defined correctly!\n");
		exit(1);
	}


	return 0;
}
