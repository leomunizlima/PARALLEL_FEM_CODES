#include "baroclinic.h"

double BAROCLINIC_epresc(double x, double y)
{
	double gamma = 1.4, p0 = 1.0, epsilon = 0.05, rho0 = 1.0, rho1 = 0.001, rho2 = 1.8;
	double rho, p1, p, Phiy, epsilon2, U4, v1, u0;
	p1 = gamma;
	u0 = sqrt( gamma );
	epsilon2 = epsilon * epsilon;
	
	if( y <= 4.0 )
			Phiy = (rho2/8.0) * y;
		else
			Phiy = rho2 * (y/8.0 - 1.0);

	if( fabs(y) <= 1e-15 || fabs(y - 8.0)  <= 1e-15 ){ 
		rho = rho0 + 0.5 * rho1 * epsilon * ( 1.0 + cos( PI * x / 20.0 ) ); //  + Phi(0), Phi(0) = Phi(8) = 0
		v1 = 0.5 * u0 * ( 1.0 + cos( PI * x / 20.0 ) );		
		p = p0 + 0.5 * epsilon * p1 * ( 1.0 + cos( PI * x / 20.0 ) );	
		U4 = p / ( gamma - 1.0 ) + 0.5 * epsilon2 * rho * ( v1 * v1 ); //state equation and U4 = rho*E. ||v||^2 = v1*v1, v2=0.
	}
	else{
		rho = rho0 + Phiy;
		v1 = 0.0;
		p = 1.0;	
		U4 = p / ( gamma - 1.0 ) + 0.5 * epsilon2 * rho * ( v1 * v1 ); 
	
	}
		
	return U4;
} 

