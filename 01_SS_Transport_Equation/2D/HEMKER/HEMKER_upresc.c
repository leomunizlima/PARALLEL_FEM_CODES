#include "../01_CommonFiles/SSTranspEquation.h"
#include "hemker.h"

inline double HEMKER_upresc(double X, double Y)
{
	double u, epsilon=1e-5;

	if (fabs(X*X + Y*Y)<=1+epsilon)
		u = 1.0;
	else
		u = 0.0;

	return u;

}






