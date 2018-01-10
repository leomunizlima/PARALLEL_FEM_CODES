#include "../01_CommonFiles/TranspEquation.h"
#include "pudim.h" 

double PUDIM_upresc(double X, double Y)
{
	double u;

	//Pudim
	if (fabs(X)<=1e-14 && Y<=0.)
		u = -sin(PI*Y);
	else
		u = 0;

	return u;

}






