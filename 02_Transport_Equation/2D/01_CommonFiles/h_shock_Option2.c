#include "TranspEquation.h"

double h_shock_Option2(double Be_x, double Be_y, double u1, double u2, double u3, double y23, double y31, double y12, double x32, double x13, double x21, double Area)
{
	double h_shock, P1, P2, P3, normBeta;

	P1 = fabs(Be_x*y23 + Be_y*x32);
	P2 = fabs(Be_x*y31 + Be_y*x13);
	P3 = fabs(Be_x*y12 + Be_y*x21);
	
	normBeta = sqrt(Be_x*Be_x + Be_y*Be_y);

	if (P1+P2+P3>=1e-12)
		h_shock = 4*Area*normBeta/(P1 + P2 + P3);
	else
		h_shock = 0;
	
	return h_shock;
}


