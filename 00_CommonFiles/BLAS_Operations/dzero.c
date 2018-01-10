#include "ourBLAS.h"

int dzero(int n, double *x)
{
	long int i, m;

      	m = n-3;
	for (i = 0; i < m; i += 4){
		x[i] = 0.;
		x[i+1] = 0.;
	        x[i+2] = 0.;
	        x[i+3] = 0.;
	}
	for ( ; i < n; ++i) 
		x[i] = 0.;
	
	return 0;
}


