#include "reaction.h" 
#include "math.h"

inline double REACTION_Font(double X, double Y, double k, double gamma, double Be_x, double Be_y)
{
	double f;

	if (X <= 0.5)
		f = X;
	else
		f = 1.0 - X;
		
	return f;

}



