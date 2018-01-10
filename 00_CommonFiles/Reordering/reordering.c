#include "reordering.h"

void reordering(ParametersType *Parameters, int *JA, int *IA, int *perm)
{
	int neq;

	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0)
		neq = Parameters->neq;
	else
		neq = Parameters->NEQ;

	if (strcasecmp(Parameters->reordering,"NOT")==0)
		return;
	else if (strcasecmp(Parameters->reordering,"Spectral")==0)
		REORDERING_SPECTRAL (neq, Parameters->nnzero, JA,  IA, perm);
	else if (strcasecmp(Parameters->reordering,"SYMRCM")==0)
		REORDERING_SYMRCM (neq, Parameters->nnzero,JA,  IA, perm);
	else{
		printf("Reordering scheme not defined!\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}
}

/*----------------------------------------------------------------------------
 * SPECTRAL reordering
 *--------------------------------------------------------------------------*/
void REORDERING_SPECTRAL (int n, int nz, int *ja, int *ia, int *p)
{
	int i;
	
	double *fvector = calloc (n ,sizeof(double));
	int    *list    = calloc (n ,sizeof(int));
	int    *info    = calloc (10,sizeof(int));
	int   lirn = nz;

		
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nz; i++) {
		ja[i] += 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] += 1;
	}

	mc73_fiedler(&n,&lirn,ja,ia,list,fvector,info,NULL);
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nz; i++) {
		ja[i] -= 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] -= 1;
	}
	
	ARRAY2* R = malloc(n*sizeof(ARRAY2));
	for (i = 0; i < n; ++i)
	{
		R[i].arr1 = fvector[i];
		R[i].arr2 = i;		
	}	
	
	qsort (R,n,sizeof(ARRAY2),COMPARE_eig); 
	
	for (i = 0; i < n; ++i) 
		p[i] = R[i].arr2; 

	free(R);
	free(info);
	free(list);
	free(fvector);
	
	return;
}

int COMPARE_eig (const void * a, const void * b)
{
	if (((ARRAY2*)a)->arr1 > ((ARRAY2*)b)->arr1) return  1;
	if (((ARRAY2*)a)->arr1 < ((ARRAY2*)b)->arr1) return -1;
	return 0;
}

int COMPARE_array (const void * a, const void * b)
{
	if (((ARRAY*)a)->arr3 <  ((ARRAY*)b)->arr3) return -1;
	if (((ARRAY*)a)->arr3 >  ((ARRAY*)b)->arr3) return  1;
	if (((ARRAY*)a)->arr3 == ((ARRAY*)b)->arr3)
	{
		if (((ARRAY*)a)->arr2 < ((ARRAY*)b)->arr2) return -1;
		if (((ARRAY*)a)->arr2 > ((ARRAY*)b)->arr2) return  1;
	}
	return 0;
}

/*----------------------------------------------------------------------------
 * Perform the colunm permutation
 *--------------------------------------------------------------------------*/
void MATRIX_COL_permutation (int n, int nz, int *JA, int *IA, int *p, int *pT)
{
	int i, j, k;

	ARRAY* a = calloc (nz,sizeof(ARRAY));
	int*   q = calloc (n ,sizeof(int));
	
	for (i = 0; i < n; ++i) 
		q[p[i]] = i; 

	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = IA[i]; j <= IA[i+1] - 1; ++j)
		{
			a[k].arr1 = pT[j];
			a[k].arr2 = q[JA[j]];
			a[k].arr3 = i;
				k = k + 1;
		}
		IA[i+1] = k;    
	}

	qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	
	for (i = 0; i < nz; ++i)
	{
		pT[i] = a[i].arr1;
		JA[i] = a[i].arr2;
	}

	free(a);
	free(q);
}

/*----------------------------------------------------------------------------
 * Perform the colunm permutation
 *--------------------------------------------------------------------------*/
void MATRIX_ROW_permutation (int n, int nz, int *JA, int *IA, int *p, int *pT)
{
	int i, j, k;
	
	int* 	auxpT = malloc( nz  *sizeof(int));
	int*    auxJA = malloc( nz  *sizeof(int));
	int*    auxIA = malloc((n+1)*sizeof(int));
  
	auxIA[0] = 0;
	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = IA[p[i]]; j <= IA[p[i]+1] - 1; ++j)
		{
			auxpT[k] = pT[j];
			auxJA[k] = JA[j];
			      k  = k + 1;
		}
		auxIA[i+1] = k;    
	}

	memcpy(&pT[0],&auxpT[0],nz*sizeof(int));
	memcpy(&JA[0],&auxJA[0],nz*sizeof(int));
	memcpy(&IA[0],&auxIA[0],(n+1)*sizeof(int));

	free(auxpT);
	free(auxJA);
	free(auxIA);
}

void MATRIX_SPARSE_ROW_permutation (int n, int nz, int *JA, int *IA, int *aux, int *pT)
{
	int i, j, k;
	
	int* 	auxpT = calloc( nz  ,sizeof(int));
	int*    auxJA = calloc( nz  ,sizeof(int));
	int*    auxIA = calloc((n+1),sizeof(int));

	ARRAY* Temp = calloc(n, sizeof(ARRAY));  	

	for (i = 0; i < n; i++)
	{
		Temp[i].arr3 = aux[i];
		Temp[i].arr2 = i;
	}

	qsort (Temp,n,sizeof(ARRAY),COMPARE_array); 

	auxIA[0] = 0;
	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = IA[Temp[i].arr2]; j <= IA[Temp[i].arr2+1] - 1; ++j)
		{
			auxpT[k] = pT[j];
			auxJA[k] = JA[j];
			      k  = k + 1;
		}
		auxIA[i+1] = k;    
	}


	for (i=0; i<n; i++)
		aux[i] = Temp[i].arr3;

	free(Temp);

	memcpy(&pT[0],&auxpT[0],nz*sizeof(int));
	memcpy(&JA[0],&auxJA[0],nz*sizeof(int));
	memcpy(&IA[0],&auxIA[0],(n+1)*sizeof(int));

	free(auxpT);
	free(auxJA);
	free(auxIA);
}


int MATRIX_bandwidth (int n, int *JA, int *IA)
{
	int i;
	int bandl, bandr;
	int bandwidth=0;

	bandl = 0;
	bandr = 0;

	for (i = 0; i < n ; i++)
	{
		if (fabs(i - JA[IA[i]]) > bandl)
			bandl= i - JA[IA[i]];
		if (fabs(JA[IA[i+1]-1]-i) > bandr);	
			bandr = JA[IA[i+1]-1]-i; 
	}
	
	if (bandl>bandr)
		bandwidth = bandl;
	else
		bandwidth = bandr;

	return bandwidth;
	
}


