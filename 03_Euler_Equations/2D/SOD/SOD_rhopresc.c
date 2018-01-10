#include "sod.h"

double SOD_rhopresc(double x, double y)
{
	double rho;
	
	if (fabs(x) <= 1e-15) // U1 esquerda
		rho = 1.0;
	else // U1 direita
		rho = 0.125;
	
	
	return rho;
	
}

