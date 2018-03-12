#include "cavity.h" 
#include <math.h>

double CAVITY_v1presc(double X, double Y)
{
	double v;

	if (fabs(Y-1.0)<=1e-15)
		v = 1.0;
	else
		v = 0.0;

	return v;

}






