#include "refletido.h"

double REFLETIDO_v2presc(double x, double y){
	double v2;
	
	if (fabs(x)<=1e-15 || fabs(y)<=1e-15) // rho*v2 - esquerda ou rho*v2 - inferior
		v2 = 0.0;
	else  // rho*v2 - topo
		v2 = -0.860744;
	
	
	return v2;
} 



