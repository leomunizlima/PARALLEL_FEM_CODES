#include "TranspEquation.h"

double h_shock_Option1(double Be_x, double Be_y, double u1, double u2, double u3, double y23, double y31, double y12, double x32, double x13, double x21, double Area)
{
	double h_shock, P1, P2, P3, aux1, aux2, normGraduArea;

	P1 = fabs(u1*y23*y23 + u2*y31*y23 + u3*y12*y23 + u1*x32*x32 + u2*x13*x32 + u3*x21*x32);
	P2 = fabs(u1*y23*y31 + u2*y31*y31 + u3*y12*y31 + u1*x32*x13 + u2*x13*x13 + u3*x21*x13);
	P3 = fabs(u1*y23*y12 + u2*y31*y12 + u3*y12*y12 + u1*x32*x21 + u2*x13*x21 + u3*x21*x21);
	aux1 = y23*u1 + y31*u2 + y12*u3;
   	aux2 = x32*u1 + x13*u2 + x21*u3;
	normGraduArea = sqrt(aux1*aux1 + aux2*aux2);
	if (P1+P2+P3>= 1e-12)
		h_shock = 4*Area*normGraduArea/(P1 + P2 + P3);
	else
		h_shock = 0;
	
	return h_shock;
}


