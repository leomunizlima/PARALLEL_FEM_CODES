#include "ourBLAS.h"

int dcopy(int n, double *x, double *y)
{
	long int i, m;

      	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i] = x[i];
		y[i+1] = x[i+1];
	        y[i+2] = x[i+2];
	        y[i+3] = x[i+3];
	}
	for ( ; i < n; ++i) /* clean-up loop */
		y[i] = x[i];
	
	return 0;
}


