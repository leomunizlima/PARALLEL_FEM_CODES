#include "cone2.h"

void CONE2_Velocity(double X, double Y, double Beta[2])
{
	double Be_x, Be_y;

	Be_x = -Y;
	Be_y = X;
	Beta[0] = Be_x;
	Beta[1] = Be_y;	

}


