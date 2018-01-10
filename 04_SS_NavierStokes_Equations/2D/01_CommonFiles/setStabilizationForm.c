#include "SSNavierStokesEquations.h"

int setStabilizationForm(ParametersType *Parameters,FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	
	if (strcasecmp(Parameters->StabilizationForm,"SUPG-PSPG")==0){
			FemOtherFunctions->Build = Build_K_F_SUPG_PSPG;
	}		
	else if (strcasecmp(Parameters->StabilizationForm,"VMS-DCDD")==0){
			FemOtherFunctions->Build = Build_K_F_VMS_DCDD;
	}		
	else {
		printf("Stabilizantion form not defined!\n");
		exit(1);
	}

	return 0;
}

