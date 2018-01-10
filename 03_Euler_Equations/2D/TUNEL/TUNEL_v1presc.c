#include "tunel.h"

double TUNEL_v1presc(double x, double y)
{
	double v1;

	if(fabs(x)<=1e-15) // U2 = rho*v1 = 1.4 * 3.0 - esquerda
		v1 = 4.2;
	else // U2 = rho*v1 = 1.4 * 0.0 - esquerda do step
		v1 = 0.0;
	
	return v1;

}
 
