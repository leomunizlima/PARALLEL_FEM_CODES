#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int SGS_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{

	SGS_precond_EBE (Parameters, MatrixData, FemStructs, &F[Parameters->neq_bef], &F[Parameters->neq_bef]);

	return 0;

}

