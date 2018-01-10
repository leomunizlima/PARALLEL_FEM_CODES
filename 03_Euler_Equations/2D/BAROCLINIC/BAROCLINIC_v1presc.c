#include "baroclinic.h"

double BAROCLINIC_v1presc(double x, double y)
{
	double v1, u0, rho, U2, Phiy; 
	double gamma = 1.4, rho0 = 1.0, rho1 = 0.001, rho2 = 1.8, epsilon = 0.05;
	u0 = sqrt( gamma );

	if( y <= 4.0 )
			Phiy = (rho2/8.0) * y;
		else
			Phiy = rho2 * (y/8.0 - 1.0);

	if( fabs(y) <= 1e-15 || fabs(y - 8.0)  <= 1e-15 ){
 		rho = rho0 + 0.5 * rho1 * epsilon * ( 1.0 + cos( PI * x / 20.0 ) ); //  + Phi(0), Phi(0) = Phi(8) = 0
		v1 = 0.5 * u0 * ( 1.0 + cos( PI * x / 20.0 ) );
	}
	else{
 		rho = rho0 + Phiy;
		v1 = 0.0;
	}
	U2 = rho * v1;	

	return U2;

}
 
