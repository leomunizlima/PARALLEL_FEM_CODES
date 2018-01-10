#include "ourBLAS.h"

int dtrsvUP (int n, double **A, double *b, double *x)
{
	int i,j;
	double sum;

	for (i=n-1;i>=0;i--)
	{
		sum = 0;
		for (j=i+1;j<=n-1;j++)
			sum += A[j][i]*x[j];
		x[i] = (1./A[i][i])*(b[i]-sum);
	}
	
	return 0;
}



