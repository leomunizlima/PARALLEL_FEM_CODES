#include "TranspEquation.h"

void csr_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double ke[3][3])
{
	double *K, *K_bef, *K_aft;
	int nnzero, nnzero_bef, nnzero_aft, *CSR_by_Element, *CSR_by_Element_bef, *CSR_by_Element_aft;
	
	nnzero = Parameters->nnzero;
	nnzero_bef = Parameters->nnzero_bef;
	nnzero_aft = Parameters->nnzero_aft;
	K = MatrixData->AA;	
	K_bef = MatrixData->AA_bef;	
	K_aft = MatrixData->AA_aft;	
	CSR_by_Element = MatrixData->Scheme_by_Element[E];
	CSR_by_Element_bef = MatrixData->Scheme_by_Element_bef[E];
	CSR_by_Element_aft = MatrixData->Scheme_by_Element_aft[E];
	
	K[CSR_by_Element[0]] += ke[0][0];
	K[CSR_by_Element[1]] += ke[0][1];
	K[CSR_by_Element[2]] += ke[0][2];
	K[CSR_by_Element[3]] += ke[1][0];
	K[CSR_by_Element[4]] += ke[1][1];
	K[CSR_by_Element[5]] += ke[1][2];
	K[CSR_by_Element[6]] += ke[2][0];
	K[CSR_by_Element[7]] += ke[2][1];
	K[CSR_by_Element[8]] += ke[2][2];
 
	K[nnzero] = 0;
	
	K_bef[CSR_by_Element_bef[0]] += ke[0][0];
	K_bef[CSR_by_Element_bef[1]] += ke[0][1];
	K_bef[CSR_by_Element_bef[2]] += ke[0][2];
	K_bef[CSR_by_Element_bef[3]] += ke[1][0];
	K_bef[CSR_by_Element_bef[4]] += ke[1][1];
	K_bef[CSR_by_Element_bef[5]] += ke[1][2];
	K_bef[CSR_by_Element_bef[6]] += ke[2][0];
	K_bef[CSR_by_Element_bef[7]] += ke[2][1];
	K_bef[CSR_by_Element_bef[8]] += ke[2][2];

	K_bef[nnzero_bef] = 0;
	
	K_aft[CSR_by_Element_aft[0]] += ke[0][0];
	K_aft[CSR_by_Element_aft[1]] += ke[0][1];
	K_aft[CSR_by_Element_aft[2]] += ke[0][2];
	K_aft[CSR_by_Element_aft[3]] += ke[1][0];
	K_aft[CSR_by_Element_aft[4]] += ke[1][1];
	K_aft[CSR_by_Element_aft[5]] += ke[1][2];
	K_aft[CSR_by_Element_aft[6]] += ke[2][0];
	K_aft[CSR_by_Element_aft[7]] += ke[2][1];
	K_aft[CSR_by_Element_aft[8]] += ke[2][2];

	K_aft[nnzero_aft] = 0;

}




