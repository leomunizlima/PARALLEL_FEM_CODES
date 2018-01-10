#include "cartola.h" 

void CARTOLA_Velocity(double X, double Y, double Beta[2])
{
	double Be_x, Be_y;
	
	double xo=0.5, yo=0.5, ro=0.25;
	
	if ((X-xo)*(X-xo)-(Y-yo)*(Y-yo)<=ro*ro && (X-xo)*(X-xo)-(Y-yo)*(Y-yo)>=0){
		Be_x = -(2*Y-1)*(ro*ro -(X-xo)*(X-xo)-(Y-yo)*(Y-yo));
		Be_y = (2*X-1)*(ro*ro -(X-xo)*(X-xo)-(Y-yo)*(Y-yo));
	}
	else{
		Be_x = 0;
		Be_y = 0;
	}

	Beta[0] = Be_x;
	Beta[1] = Be_y;	

}


