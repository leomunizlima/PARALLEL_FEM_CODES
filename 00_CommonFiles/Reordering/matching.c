/*----------------------------------------------------------------------------
 * MATCHING REORDERING
 *--------------------------------------------------------------------------*/
#include "../COMMON_FILES/protos.h"

/*----------------------------------------------------------------------------
 * MATCHING reordering
 *--------------------------------------------------------------------------*/
void mc64ad_(int* JOB,int* N,int* NE,int* IP,int* IRN,double* AA,int* NUM,int* CPERM,int* LIW,int* IW,int* LDW,double* DW,int* ICNTL,int* INFO);

void REORDERING_MATCHING (MAT* A, double* b, int** Q, double** dw, int* msgs)
{
	double time;
	/*---START TIME---------------> */ time =  get_TIME(); 
	
	int i,j,k;
	int JOB = 5;
	int N = A->n;
	int NE = A->nz;
	int NUM;
	int MAXN  = N; 
	int LIW   = 100*N;
	int LDW   = 100*N;
	
	int* IP = A->IA;
	int* IRN = A->JA;
	int* CPERM = calloc (N,sizeof(int));
	int* IW = calloc (LIW,sizeof(int));
	int* ICNTL = calloc (10,sizeof(int));
	int* INFO = calloc (10,sizeof(int));
	
	double* AA = A->AA;
	double* DW = calloc (LDW,sizeof(double));
	
	ICNTL[0] = 0;
        ICNTL[1] = -1;
        ICNTL[2] = -1;
        for (i = 3;i < 10; ++i)
		ICNTL[i] = 0;

	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < NE; i++) {
		IRN[i] += 1;
	}
	for (i = 0; i < N + 1; i++) {
		IP[i] += 1;
	}

	mc64ad_(&JOB,&N,&NE,IP,IRN,AA,&NUM,CPERM,&LIW,IW,&LDW,DW,ICNTL,INFO);
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < NE; i++) {
		IRN[i] -= 1;
	}
	for (i = 0; i < N + 1; i++) {
		IP[i] -= 1;
	}
	
	/* -------------------------------------------------------------------- */    
	/* Appling scaling and matching permutations                            */
	/* -------------------------------------------------------------------- */
	if (msgs[1])
	{
		for (j = 0; j < N; ++j)
		{
			for (k = A->IA[j]; k <= A->IA[j+1] - 1; ++k)
			{
				A->AA[k] = A->AA[k] * exp(DW[A->JA[k]] + DW[N+j]);
			}
			b[j] = b[j] * exp(DW[N+j]);
		}
	}
	
	if (msgs[2])
	{
		for (i = 0; i < N; ++i)
			CPERM[i]--;	
	
		MATRIX_ROW_permutation (A,CPERM);
		
			/*---PERMUTE COEF. VECTOR------ */
		double* bp = calloc(A->n,sizeof(double));		
		for (i = 0; i < A->n; ++i) 
			bp[i] = b[CPERM[i]];
		
		memcpy(&b[0],&bp[0],A->n*sizeof(double));
		free(bp);
		
	}
	
	(*Q)  = CPERM;
	(*dw) = DW;
	
	free(INFO);
	free(ICNTL);
	free(IW);
	
	/*---FINAL TIME---------------> */ time = (get_TIME() - time)/100.0;
	
	if (msgs[0]) fprintf(stderr,"\n  [ MATCHING AND SCALING ]\n");
	if (msgs[0]) fprintf(stderr,"  - Elapsed time: %.6f sec\n", time);
	
	return;
}
