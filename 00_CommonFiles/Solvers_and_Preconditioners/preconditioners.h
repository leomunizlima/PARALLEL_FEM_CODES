#ifdef SSTranspEquation2D
	#include "../../01_SS_Transport_Equation/2D/01_CommonFiles/SSTranspEquation.h"
	int LU_precond_EBE (ParametersType *, MatrixDataType *, FemStructsType *,double *, double *);
	int LU_precond_EBE_setup (ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
	int SGS_precond_EBE (ParametersType *, MatrixDataType *, FemStructsType *,double *, double *);
	int SGS_precondR_EBE (ParametersType *, MatrixDataType *, FemStructsType *,double *, double *);
	int SGS_precond_EBE_setup (ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
#endif
#ifdef TranspEquation2D
	#include "../../02_Transport_Equation/2D/01_CommonFiles/TranspEquation.h"
	int LU_precond_EBE (ParametersType *, MatrixDataType *, FemStructsType *,double *, double *);
	int LU_precond_EBE_setup (ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
	int SGS_precond_EBE (ParametersType *, MatrixDataType *, FemStructsType *,double *, double *);
	int SGS_precondR_EBE (ParametersType *, MatrixDataType *, FemStructsType *,double *, double *);
	int SGS_precond_EBE_setup (ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
#endif
#ifdef EulerEquations2D
	#include "../../03_Euler_Equations/2D/01_CommonFiles/EulerEquations.h"
	int BlockDiag_precond (ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int BlockDiag_precond_EBE_setup(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
#endif
#ifdef SSNavierStokesEquations2D
	#include "../../04_SS_NavierStokes_Equations/2D/01_CommonFiles/SSNavierStokesEquations.h"
	int BlockDiagDOF3_precond (ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int BlockDiagDOF3_precond_EBE_setup(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
#endif

int NO_precond   (ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
int Diag_precond   (ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
int ILUp_precond   (ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
int SPIKE_precond   (ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
int NO_precond_setup(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
int Diag_precond_EBE_setup(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
int Diag_precond_CSR_setup(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
int ILUp_precond_setup(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
int SPIKE_precond_setup(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);


