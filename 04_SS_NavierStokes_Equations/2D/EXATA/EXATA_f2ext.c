#include "exata.h" 
#include <math.h>

double EXATA_f2ext(double x, double y)
{
	double f;

	f = -2*y - 1.2*pow((1 - y),2)*pow(y,2) + 8*pow(x,2)*pow((1 - 3*x + 2*pow(x,2)),2)*pow((-1 + y),3)*pow(y,3)*(-1 + 2*y) + pow(x,2)*(-1.2 + 7.2*y - 7.2*pow(y,2)) - 4*pow((-1 + x),2)*pow(x,2)*(1 - 6*x + 6*pow(x,2))*pow((-1 + y),2)*pow(y,3)*(1 - 3*y + 2*pow(y,2)) + pow(x,3)*(0.8 - 4.8*y + 4.8*pow(y,2)) + x*(0.4 - 2.4*y + 4.8*pow(y,2) - 4.8*pow(y,3) + 2.4*pow(y,4));

//	f = 1e-1;

	return f;

}






