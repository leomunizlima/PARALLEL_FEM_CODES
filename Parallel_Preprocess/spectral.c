#include "protos.h"

/*----------------------------------------------------------------------------
 * SPECTRAL reordering
 *--------------------------------------------------------------------------*/
void mc73_fiedler(int* n, int* lirn, int* irn, int* ip, int* list, double* fvector, int* info, double* a);
int COMPARE_eig (const void * a, const void * b);

void REORDERING_SPECTRAL (int n, int nnz, int *ja, int *ia, int *p)
{
	int i;		
	double *fvector = calloc (n ,sizeof(double));;
	int    *list    = calloc (n ,sizeof(int));
	int    *info    = calloc (10,sizeof(int));
	int   lirn = nnz;
	
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] += 1;
	}

	mc73_fiedler(&n,&lirn,ja,ia,list,fvector,info,NULL);
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] -= 1;
	}
	
	ARRAY* R = malloc(n*sizeof(ARRAY));
	for (i = 0; i < n; ++i)
	{
		R[i].arr1 = fvector[i];
		R[i].arr2 = i;		
	}
	
	qsort (R,n,sizeof(ARRAY),COMPARE_eig); 
	
	for (i = 0; i < n; ++i) 
		p[i] = R[i].arr2; 

	free(R);	
	free(fvector);
	free(info);
	free(list);

	return;
}

int COMPARE_eig (const void * a, const void * b)
{
	if (((ARRAY*)a)->arr1 > ((ARRAY*)b)->arr1) return  1;
	if (((ARRAY*)a)->arr1 < ((ARRAY*)b)->arr1) return -1;
	return 0;
}


