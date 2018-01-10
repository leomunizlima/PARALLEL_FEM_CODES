#ifdef TranspEquation2D
	#include "../../02_Transport_Equation/2D/01_CommonFiles/TranspEquation.h"
#endif
#ifdef EulerEquations2D
	#include "../../03_Euler_Equations/2D/01_CommonFiles/EulerEquations.h"
#endif
#include "../BLAS_Operations/ourBLAS.h"
#include "../Allocation_Operations/allocations.h"
#include "../IO_Operations/io.h"
#include "../C_Operations/c.h"

int Predictor_Old(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_Old_BDF(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New_BDF(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New_TRBDF2(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_Old_TRBDF2(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

void setStopCriteria(ParametersType *, FemFunctionsType *);

int StopByIterations(ParametersType *, double, double, int);

int StopByNorm(ParametersType *, double, double, int);

int StopBySteadyState(ParametersType *, double *, double *, double);

int StopByTime(ParametersType *, double *, double *, double);







