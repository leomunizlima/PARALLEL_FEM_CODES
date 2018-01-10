#include "c.h"

void dmemcpy(int n, double *x, double *y)
{
	memcpy(y, x, n*sizeof(double)); //copy x to y
}


