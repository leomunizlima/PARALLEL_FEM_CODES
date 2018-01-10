#include "../01_CommonFiles/SSTranspEquation.h"
#include "bigpudim.h" 

inline double BIGPUDIM_upresc(double X, double Y)
{
	double u;

	//Pudim
	if (fabs(X)<=1e-14 && Y<=0.)
		u = -10*sin(0.1*PI*Y);
	else
		u = 0;

	return u;

}






