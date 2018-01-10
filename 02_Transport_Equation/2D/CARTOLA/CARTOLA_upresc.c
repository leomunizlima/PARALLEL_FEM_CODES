#include "../01_CommonFiles/TranspEquation.h"
#include "cartola.h" 

double CARTOLA_upresc(double X, double Y)
{
	double u;
	double A, a = 1000, ro = 0.25, xo = 0.5, yo = 0.5 ;

	A = a * ( ro*ro - (X-xo)*(X-xo) - (Y-yo)*(Y-yo) );
	u = 0.5 + atan(A)/PI;


	return u;

}






