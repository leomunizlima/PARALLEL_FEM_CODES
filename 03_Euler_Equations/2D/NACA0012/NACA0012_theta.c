#include "naca0012.h"

double NACA0012_theta(double x, double y) 
{
	double theta;

	if (fabs(0.178140*sqrt(x) - 0.075600*x - 0.210960*x*x +  0.170580*x*x*x - 0.060900*x*x*x*x -y)<=1e-10)
		theta = atan(0.089070/sqrt(x) - 0.075600 - 0.421920*x + 0.511740*x*x - 0.243600*x*x*x);
	else
		theta = atan(-0.089070/sqrt(x) + 0.075600 + 0.421920*x - 0.511740*x*x + 0.243600*x*x*x);


	return theta;
}

