#include "SSTranspEquation.h"

int setStabilizationForm(ParametersType *Parameters, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions) 
{
	if (strcasecmp(Parameters->StabilizationForm,"SUPG")==0){
	
		FemOtherFunctions->Build_K_F = Build_K_F_SUPG;

		if (strcasecmp(Parameters->ShockCapture,"CAU")==0){
		 	FemFunctions->ShockCapture = CAU_ShockCapture;
		}
		else if (strcasecmp(Parameters->ShockCapture,"YZBeta")==0){
			FemFunctions->ShockCapture = YZBeta_ShockCapture;
		}	
		else{
			printf("Shock capture not defined correctly!\n");
			MPI_Finalize();
		}
		
	}
	else if (strcasecmp(Parameters->StabilizationForm,"DD")==0){
	
		FemOtherFunctions->Build_K_F = Build_K_F_DD;
		
		if (strcasecmp(Parameters->ShockCapture,"CAUDD")==0){
			FemFunctions->ShockCapture = CAU_DD_ShockCapture;
		}
		else if (strcasecmp(Parameters->ShockCapture,"CAU")==0){
			FemFunctions->ShockCapture = CAU_ShockCapture;
		}
		else if (strcasecmp(Parameters->ShockCapture,"YZBeta")==0){
			FemFunctions->ShockCapture = YZBeta_ShockCapture;
		}	
		else{
			printf("Shock capture not defined correctly!\n");
			MPI_Finalize();
		}	
	}
	else {
		printf("Stabilizantion form not defined!\n");
		MPI_Finalize();
	}

	if (strcasecmp(Parameters->h_Shock,"2sqrtArea")==0){ 
		FemFunctions->h_shock_calculation = h_shock_2sqrtArea;
	}
	else if (strcasecmp(Parameters->h_Shock,"Option1")==0){  
		FemFunctions->h_shock_calculation = h_shock_Option1;
	}
	else if (strcasecmp(Parameters->h_Shock,"Option2")==0){ 
		FemFunctions->h_shock_calculation = h_shock_Option2;
	}
 	else{
		printf("h_shock not defined!\n");
		MPI_Finalize();
	}



	return 0;
}


