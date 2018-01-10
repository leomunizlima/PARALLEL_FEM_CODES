#include "ourBLAS.h"

int izero(int n, int *v)
{
	long int i, m;

      	m = n-3;
	for (i = 0; i < m; i += 4){
		v[i] = 0;
		v[i+1] = 0;
	        v[i+2] = 0;
	        v[i+3] = 0;
	}
	for ( ; i < n; ++i) 
		v[i] = 0;
	
	return 0;
}


