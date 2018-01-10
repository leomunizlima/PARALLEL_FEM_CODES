#include "TranspEquation.h"
#include "../../../00_CommonFiles/Time_Integration/time_integration.h"	

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, 
	    FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int (*Predictor)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

	setProblem(Parameters, FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters,FemOtherFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions, &Predictor);
	setStopCriteria(Parameters, FemFunctions);
	Predictor(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
	
	return 0;
}


