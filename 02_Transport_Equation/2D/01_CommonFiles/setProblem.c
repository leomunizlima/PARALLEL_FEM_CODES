#include "TranspEquation.h"
#include "../PUDIM/pudim.h"
#include "../CARTOLA/cartola.h"
#include "../CONE/cone.h"
#include "../CONE2/cone2.h"

int setProblem(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{		
	if (strcasecmp(Parameters->ProblemTitle,"PUDIM")==0){
		FemFunctions->Condutivity = PUDIM_Condutivity;	
		FemFunctions->Font = PUDIM_Font;
		FemFunctions->Reaction = PUDIM_Reaction;
		FemFunctions->Velocity = PUDIM_Velocity;
		FemFunctions->upresc = PUDIM_upresc;
		FemFunctions->InitialSolution = PUDIM_InitialSolution;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"CARTOLA")==0){
		FemFunctions->Condutivity = CARTOLA_Condutivity;	
		FemFunctions->Font = CARTOLA_Font;
		FemFunctions->Reaction = CARTOLA_Reaction;
		FemFunctions->Velocity = CARTOLA_Velocity;
		FemFunctions->upresc = CARTOLA_upresc;
		FemFunctions->InitialSolution = CARTOLA_InitialSolution;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"CONE")==0){
		FemFunctions->Condutivity = CONE_Condutivity;	
		FemFunctions->Font = CONE_Font;
		FemFunctions->Reaction = CONE_Reaction;
		FemFunctions->Velocity = CONE_Velocity;
		FemFunctions->upresc = CONE_upresc;
		FemFunctions->InitialSolution = CONE_InitialSolution;
	}else if (strcasecmp(Parameters->ProblemTitle,"CONE2")==0){
		FemFunctions->Condutivity = CONE2_Condutivity;	
		FemFunctions->Font = CONE2_Font;
		FemFunctions->Reaction = CONE2_Reaction;
		FemFunctions->Velocity = CONE2_Velocity;
		FemFunctions->upresc = CONE2_upresc;
		FemFunctions->InitialSolution = CONE2_InitialSolution;
	}
	else{
		printf("Problem not defined!\n");
		MPI_Finalize();
	}

	
	return 0;
}



