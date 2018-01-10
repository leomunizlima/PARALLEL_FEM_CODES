#include "cartola.h" 

double CARTOLA_Font(double X, double Y, double k, double gamma, double Be_x, double Be_y)
{
	double f;
	double A, a = 1000, ro = 0.25, xo = 0.5, yo = 0.5 ;

	A = a * ( ro*ro - (X-xo)*(X-xo) - (Y-yo)*(Y-yo) );
	f = -4*a*k/(PI*(A*A+1)) -(2*a*Be_x*(X-xo) + 2*a*Be_y*(Y-yo))/(PI*(A*A+1)) + gamma*(0.5 + atan(A)/PI);

	return f;

}



