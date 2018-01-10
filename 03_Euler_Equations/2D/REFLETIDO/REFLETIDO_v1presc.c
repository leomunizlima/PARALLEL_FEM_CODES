#include "refletido.h"

double REFLETIDO_v1presc(double x, double y){
	double v1;

	if(fabs(x)<=1e-15) // rho*v1 - esquerda
		v1 = 2.9;
	else  // rho*v1 - topo
		v1 = 4.452878;
	
	
	return v1;

}
 
