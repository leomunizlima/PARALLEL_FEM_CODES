#ifdef SSTranspEquation2D
	#include "../../01_SS_Transport_Equation/2D/01_CommonFiles/SSTranspEquation.h"
#endif
#ifdef TranspEquation2D
	#include "../../02_Transport_Equation/2D/01_CommonFiles/TranspEquation.h"
#endif
#ifdef EulerEquations2D
	#include "../../03_Euler_Equations/2D/01_CommonFiles/EulerEquations.h"
#endif
#ifdef SSNavierStokesEquations2D
	#include "../../04_SS_NavierStokes_Equations/2D/01_CommonFiles/SSNavierStokesEquations.h"
#endif

int pgmres (ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType*, double *B, double *X);


