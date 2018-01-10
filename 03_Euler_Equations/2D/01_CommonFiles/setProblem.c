#include "EulerEquations.h"
#include "../OBLIQUO/obliquo.h"
#include "../REFLETIDO/refletido.h"
#include "../TUNEL/tunel.h"
#include "../SOD/sod.h"
#include "../EXPLOSION/explosion.h"
#include "../BAROCLINIC/baroclinic.h"
#include "../NACA0012/naca0012.h"
#include "../CYLINDER/cylinder.h"

int setProblem(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->ProblemTitle,"OBLIQUO")==0){
		FemFunctions->gamma = OBLIQUO_gamma;
		FemFunctions->cv = OBLIQUO_cv;
		FemFunctions->rhopresc = OBLIQUO_rhopresc;
		FemFunctions->v1presc = OBLIQUO_v1presc;
		FemFunctions->v2presc = OBLIQUO_v2presc;
		FemFunctions->epresc = OBLIQUO_epresc;
		FemFunctions->InitialSolution = OBLIQUO_InitialSolution;
	}
	else if(strcasecmp(Parameters->ProblemTitle,"EXPLOSION")==0){
		FemFunctions->gamma = EXPLOSION_gamma;
		FemFunctions->cv = EXPLOSION_cv;
		FemFunctions->rhopresc = EXPLOSION_rhopresc;
		FemFunctions->v1presc = EXPLOSION_v1presc;
		FemFunctions->v2presc = EXPLOSION_v2presc;
		FemFunctions->epresc = EXPLOSION_epresc;
		FemFunctions->InitialSolution = EXPLOSION_InitialSolution;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"REFLETIDO")==0){
		FemFunctions->gamma = REFLETIDO_gamma;
		FemFunctions->cv = REFLETIDO_cv;
		FemFunctions->rhopresc = REFLETIDO_rhopresc;
		FemFunctions->v1presc = REFLETIDO_v1presc;
		FemFunctions->v2presc = REFLETIDO_v2presc;
		FemFunctions->epresc = REFLETIDO_epresc;
		FemFunctions->InitialSolution = REFLETIDO_InitialSolution;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"TUNEL")==0){
		FemFunctions->gamma = TUNEL_gamma;
		FemFunctions->cv = TUNEL_cv;
		FemFunctions->rhopresc = TUNEL_rhopresc;
		FemFunctions->v1presc = TUNEL_v1presc;
		FemFunctions->v2presc = TUNEL_v2presc;
		FemFunctions->epresc = TUNEL_epresc;
		FemFunctions->InitialSolution = TUNEL_InitialSolution;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"SOD")==0){
		FemFunctions->gamma = SOD_gamma;
		FemFunctions->cv = SOD_cv;
		FemFunctions->rhopresc = SOD_rhopresc;
		FemFunctions->v1presc = SOD_v1presc;
		FemFunctions->v2presc = SOD_v2presc;
		FemFunctions->epresc = SOD_epresc;
		FemFunctions->InitialSolution = SOD_InitialSolution;
	}else if (strcasecmp(Parameters->ProblemTitle,"BAROCLINIC")==0){
		FemFunctions->gamma = BAROCLINIC_gamma;
		FemFunctions->cv = BAROCLINIC_cv;
		FemFunctions->rhopresc = BAROCLINIC_rhopresc;
		FemFunctions->v1presc = BAROCLINIC_v1presc;
		FemFunctions->v2presc = BAROCLINIC_v2presc;
		FemFunctions->epresc = BAROCLINIC_epresc;
		FemFunctions->InitialSolution = BAROCLINIC_InitialSolution;
	}else if (strcasecmp(Parameters->ProblemTitle,"NACA0012")==0){
		FemFunctions->gamma = NACA0012_gamma;
		FemFunctions->cv = NACA0012_cv;
		FemFunctions->rhopresc = NACA0012_rhopresc;
		FemFunctions->v1presc = NACA0012_v1presc;
		FemFunctions->v2presc = NACA0012_v2presc;
		FemFunctions->epresc = NACA0012_epresc;
		FemFunctions->InitialSolution = NACA0012_InitialSolution;
	}else if (strcasecmp(Parameters->ProblemTitle,"CYLINDER")==0){
		FemFunctions->gamma = CYLINDER_gamma;
		FemFunctions->cv = CYLINDER_cv;
		FemFunctions->rhopresc = CYLINDER_rhopresc;
		FemFunctions->v1presc = CYLINDER_v1presc;
		FemFunctions->v2presc = CYLINDER_v2presc;
		FemFunctions->epresc = CYLINDER_epresc;
		FemFunctions->InitialSolution = CYLINDER_InitialSolution;
	}
	else{
		printf("Problem is not defined correctly!\n");
		exit(1);

	}
	return 0;
}


