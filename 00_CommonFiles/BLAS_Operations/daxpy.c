#include "ourBLAS.h"

int daxpy(int n, double a, double *x, double *y)
{
	long int i, m;
	register double sa;

	sa = a;
      	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i] += sa * x[i];
		y[i+1] += sa * x[i+1];
	        y[i+2] += sa * x[i+2];
	        y[i+3] += sa * x[i+3];
	}
	for ( ; i < n; ++i) /* clean-up loop */
		y[i] += sa * x[i];
	
	return 0;
}


