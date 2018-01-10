#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/solvers.h"

int setSolver(ParametersType *Parameters, FemOtherFunctionsType *FemOtherFunctions)
{
	if (strcasecmp(Parameters->Solver,"GMRES")==0){
		FemOtherFunctions->solver = pgmres;	
	}
	else{
		printf("Solver not defined correctly!\n");
		exit(1);
	}

	return 0;
}


