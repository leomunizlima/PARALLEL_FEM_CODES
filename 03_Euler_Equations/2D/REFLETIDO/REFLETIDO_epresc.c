#include "refletido.h"

double REFLETIDO_epresc(double x, double y){
	double e;
	
	if(fabs(x)<=1e-15) //rho*e - esquerda
		e = 5.990715;
	else  //rho*e - topo
		e = 9.870181681;
	

	return e;
}

 
