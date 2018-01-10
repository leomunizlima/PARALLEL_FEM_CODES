#include "refletido.h"

double REFLETIDO_rhopresc(double x, double y){
	double rho;
	
	if(fabs(x) <= 1e-15) // esquerda
		rho = 1.0;
	else // topo
		rho = 1.7;
	
	
	return rho;
	
}
