#include "TranspEquation.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/preconditioners.h"

int setPreconditioner(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Preconditioner,"NOT")==0){
		FemFunctions->precond = NO_precond;
		FemFunctions->precond_setup = NO_precond_setup;
	}
	else if (strcasecmp(Parameters->Preconditioner,"Diag")==0){
		FemFunctions->precond = Diag_precond;
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
			FemFunctions->precond_setup = Diag_precond_EBE_setup;
				
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
//			FemFunctions->precond_setup = Diag_precond_EDE_setup;

		}	 
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
			FemFunctions->precond_setup = Diag_precond_CSR_setup;
		}
	
	}
	else if (strncmp(Parameters->Preconditioner,"ILU",3)==0){
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
                	FemFunctions->precond = ILUp_precond;
			FemFunctions->precond_setup = ILUp_precond_setup;
		}
		else{
			printf("Preconditioner definied only to CSR scheme\n");
			exit(1);
		}
	}	
	else if (strncmp(Parameters->Preconditioner,"SPIKE",5)==0){
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
                	FemFunctions->precond = SPIKE_precond;
			FemFunctions->precond_setup = SPIKE_precond_setup;
		}
		else{
			printf("Preconditioner definied only to CSR scheme\n");
			exit(1);
		}
	}	
	else {
		printf("Preconditioner not defined!\n");
		exit(1);

	}

	return 0;
}


