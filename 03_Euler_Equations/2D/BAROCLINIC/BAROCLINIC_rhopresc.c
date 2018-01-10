#include "baroclinic.h"

double BAROCLINIC_rhopresc(double x, double y)
{
	double rho0 = 1.0, rho1 = 0.001, rho2 = 1.8, epsilon = 0.05;
	double rho, Phiy;
	
	if( y <= 4.0 )
			Phiy = (rho2/8.0) * y;
		else
			Phiy = rho2 * (y/8.0 - 1.0);

	if( fabs(y) <= 1e-15 || fabs(y - 8.0)  <= 1e-15 ) 
		rho = rho0 + 0.5 * rho1 * epsilon * ( 1.0 + cos( PI * x / 20.0 ) ); //  + Phi(0), Phi(0) = Phi(8) = 0
	else
		rho = rho0 + Phiy;
		
	return rho;
}


