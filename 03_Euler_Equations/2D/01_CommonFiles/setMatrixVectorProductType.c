#include "EulerEquations.h"
#include "../../../00_CommonFiles/MatrixVector_Operations/matvec.h"

int setMatrixVectorProductType(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3)==0){
		FemFunctions->assembly = ebe_assembly;
		FemFunctions->mv = ebemvNDOF4;
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		FemFunctions->assembly = csr_assembly;
		FemFunctions->mv = csrmv;
	}
	else{
		printf("\nMatrix vector product scheme is not defined correctly!\n\n");
		exit(1);
	}
	return 0;
}


