#include "sod.h"

double SOD_epresc(double x, double y)
{
	double e;
	
	if(fabs(x)<=1e-15) // U4 = rho * e = 1.0 * 2.5 - esquerda
		e = 2.5;
	else // U4 = rho * e = 0.125 * 2.0 - direita
		e = 0.25;

	return e;
} 


