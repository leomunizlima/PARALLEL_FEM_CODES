#include "../01_CommonFiles/SSTranspEquation.h"
#include "reaction.h" 
#include "math.h"

inline double REACTION_upresc(double X, double Y)
{
	double u;

	if (fabs(X)<1e-14 || fabs(Y-1)<1e-14)
		u = 1.0;
	else
		u = 0;

	return u;

}






