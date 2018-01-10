#include "TranspEquation.h"
#include "../../../00_CommonFiles/MatrixVector_Operations/matvec.h"

int setMatrixVectorProductType(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
		FemFunctions->assembly = ebe_assembly;
		FemFunctions->mv = ebemv;
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		FemFunctions->assembly = csr_assembly;
		FemFunctions->mv = csrmv;
	}
	else{
		printf("\nMatrix vector product scheme is not defined!\n\n");
		exit(1);
	}
	return 0;
}


