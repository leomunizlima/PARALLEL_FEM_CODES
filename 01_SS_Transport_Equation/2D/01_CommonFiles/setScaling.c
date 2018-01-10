#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/scaling.h"

int setScaling(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Scaling,"NOT")==0){
		FemFunctions->unscaling = NO_unscaling;
		FemFunctions->scaling = NO_scaling;
	}
	else if (strcasecmp(Parameters->Scaling,"Left")==0){
		FemFunctions->unscaling = NO_unscaling;
		if (strcasecmp(Parameters->MatrixVectorProductScheme, "EBE")==0) 
			FemFunctions->scaling = Left_scaling_EBE;
		else		
			FemFunctions->scaling = Left_scaling_CSR;

	}
	else if (strcasecmp(Parameters->Scaling,"LeftRight")==0){
		FemFunctions->unscaling = Left_unscaling;
		if (strcasecmp(Parameters->MatrixVectorProductScheme, "EBE")==0) 
			FemFunctions->scaling = LeftRight_scaling_EBE;
		else		
			FemFunctions->scaling = LeftRight_scaling_CSR;

	}
	else{
		printf("Scaling not defined correctly!\n");
		exit(1);
	}


	return 0;
}
