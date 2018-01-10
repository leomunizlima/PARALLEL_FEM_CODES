#include "EulerEquations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

void print_rho_Residue(FILE *outFile, int neqrho, double t,double *R, double *R_rho, int *eqrho)
{
	int I;
	double norm_R_rho;

	for(I = 0; I < neqrho; I++)
		R_rho[I] = R[eqrho[I]];
	
	norm_R_rho = sqrt(ddot(neqrho, R_rho, R_rho));
	
	fprintf(outFile,"%lf\t%.14e\n",t,norm_R_rho);
}


