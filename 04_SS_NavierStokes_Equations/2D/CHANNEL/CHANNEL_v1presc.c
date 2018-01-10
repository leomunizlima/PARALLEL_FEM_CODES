#include "channel.h" 
#include <math.h>

double CHANNEL_v1presc(double X, double Y)
{
	double v;

	if (fabs(X)<=1e-15)
		v = 1.0;
	else
		v = 0.0;

	return v;

}






