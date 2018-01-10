#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/solvers.h"

int setSolver(ParametersType *Parameters, FemOtherFunctionsType *FemOtherFunctions) 
{
	if (strcasecmp(Parameters->Solver,"GMRES")==0){
		FemOtherFunctions->solver = pgmres;	
	}
	else if (strcasecmp(Parameters->Solver,"CG")==0){
//		FemOtherFunctions->solver = pcg;	
	}
	else{
		printf("Solver not defined!\n");
		exit(1);
	}

	return 0;
}



