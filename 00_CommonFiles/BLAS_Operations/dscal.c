#include "ourBLAS.h"

int dscal(int n, double a, double *x)
{
	long int i, m;
	register double sa;

	sa = a;
      	m = n-3;
	for (i = 0; i < m; i += 4){
		x[i] = sa * x[i];
		x[i+1] = sa * x[i+1];
	        x[i+2] = sa * x[i+2];
	        x[i+3] = sa * x[i+3];
	}
	for ( ; i < n; ++i) /* clean-up loop */
		x[i] = sa * x[i];
	
	return 0;
}
