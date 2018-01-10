#include "EulerEquations.h"

#define size NNOEL*NDOF

void csr_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*Me)[size])
{
	int I, J, K;
	double *M, *M_bef, *M_aft;
	int nnzero, nnzero_bef, nnzero_aft, *CSR_by_Element, *CSR_by_Element_bef, *CSR_by_Element_aft;
	
	nnzero = Parameters->nnzero;
	nnzero_bef = Parameters->nnzero_bef;
	nnzero_aft = Parameters->nnzero_aft;
	M = MatrixData->AA;	
	M_bef = MatrixData->AA_bef;	
	M_aft = MatrixData->AA_aft;	
	CSR_by_Element = MatrixData->Scheme_by_Element[E];
	CSR_by_Element_bef = MatrixData->Scheme_by_Element_bef[E];
	CSR_by_Element_aft = MatrixData->Scheme_by_Element_aft[E];
		
	M = MatrixData->AA;
	CSR_by_Element = MatrixData->Scheme_by_Element[E];

	for (I=0, K=0; I<size; I++){
		for (J=0; J<size; J++){
			M[CSR_by_Element[K]] += Me[I][J];
			M_bef[CSR_by_Element_bef[K]] += Me[I][J];
			M_aft[CSR_by_Element_aft[K]] += Me[I][J];
			K++;
		}
	}

	M[nnzero] = 0;
	M_bef[nnzero_bef] = 0;
	M_aft[nnzero_aft] = 0;
}





