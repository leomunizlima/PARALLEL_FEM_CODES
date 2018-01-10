#include "cone.h"

void CONE_Velocity(double X, double Y, double Beta[2])
{
	double Be_x, Be_y;

	Be_x = -(Y - 5.0);
	Be_y = X - 5.0;
	Beta[0] = Be_x;
	Beta[1] = Be_y;	

}


