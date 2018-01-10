#include "cylinder.h"

double CYLINDER_theta(double x, double y)
{
	double theta;

	if (fabs(y)>=1e-8)
		theta = atan(-x/y); 
	else{
		printf("Math Error: division by zero!\n");
		exit(1);	
	}

	return theta;
}



