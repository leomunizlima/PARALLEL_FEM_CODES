/*----------------------------------------------------------------------------
 * PARDISO CALL FUNCTION
 *--------------------------------------------------------------------------*/
#include "pardiso.h"

/*----------------------------------------------------------------------------
 * Pardiso call function
 *--------------------------------------------------------------------------*/
void PARDISO_bottom_tips (MAT* A, double* B, double* V, int nrhs)
{
	int      n = A->n;
	int    *ia = A->IA;
	int    *ja = A->JA;
	double  *a = A->AA;
	int    nnz = A->nz;
	int  mtype = 11;        /* Real unsymmetric matrix */

	void    *pt[64];

	/* Pardiso control parameters. */
	int      iparm[64];
	double   dparm[64];
	int      solver;
	int      maxfct, mnum, phase, error, msglvl;

	/* Number of processors. */
	int      num_procs;

	/* Auxiliary variables. */
	char    *var;
	int      i, j, k;

	double   ddum;              /* Double dummy */
	int      idum;              /* Integer dummy. */
	
	error  = 0;
	solver = 0; /* use sparse direct solver */
	
	pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

	/* Numbers of processors, value of OMP_NUM_THREADS */
	var = getenv("OMP_NUM_THREADS");
	if(var != NULL)
		sscanf( var, "%d", &num_procs );
	else {
		printf("Set environment OMP_NUM_THREADS to 1\n");
		exit(1);
	}
	
	iparm[2]  = num_procs;

	iparm[10] = 1;      /* no scaling  */
	iparm[12] = 1;      /* no matching */	
	maxfct = 1;         /* Maximum number of numerical factorizations.  */
	mnum   = 1;         /* Which factorization to use. */	
	msglvl = 0;         /* Print statistical information  */
	error  = 0;         /* Initialize error flag */
	
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	int* perm = calloc (n,sizeof(int));
	for (i = n - 1; i >= n - nrhs; --i) 
		perm[i] = 1;		
	iparm[30] = 1;
	
	double* aux_V = calloc(n*nrhs,sizeof(double));
	if (aux_V == NULL)	{ 
		printf ("Memory alloc failed. Exiting.. [PARDISO_bottom_tips]\n"); 
		exit(1);		
	}
	
	/* -------------------------------------------------------------------- */    
	/* ..  Reordering and Symbolic Factorization.  This step also allocates */
	/*     all memory that is necessary for the factorization.              */
	/* -------------------------------------------------------------------- */
	phase = 11; 
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,perm,&nrhs,iparm,&msglvl,&ddum,&ddum,&error,dparm);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	
	/* -------------------------------------------------------------------- */    
	/* ..  Numerical factorization.                                         */
	/* -------------------------------------------------------------------- */    
	phase = 22;
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,perm,&nrhs,iparm,&msglvl,&ddum,&ddum,&error,dparm);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	
	/* -------------------------------------------------------------------- */    
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */ 
	phase = 33;
	iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,perm,&nrhs,iparm,&msglvl,B,aux_V,&error,dparm);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	
	/* -------------------------------------------------------------------- */    
	/*  Arrange solution                                                    */
	/* -------------------------------------------------------------------- */  
	k = 0;
	for (i = n - nrhs ; i < n; i++) {
		for (j = 0; j < nrhs; j++) {
			V[k++] = aux_V[i + j*n];
		}
	}
	free(aux_V);
	free(perm);
	
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] -= 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}

	/* -------------------------------------------------------------------- */    
	/* ..  Termination and release of memory.                               */
	/* -------------------------------------------------------------------- */ 
	phase = -1;
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,&ddum,ia,ja,&idum,&nrhs,iparm,&msglvl,&ddum,&ddum,&error,dparm);
}

void PARDISO_top_tips (MAT* A, double* C, double* W, int nrhs)
{
	int      n = A->n;
	int    *ia = A->IA;
	int    *ja = A->JA;
	double  *a = A->AA;
	int    nnz = A->nz;
	int  mtype = 11;        /* Real unsymmetric matrix */

	void    *pt[64];

	/* Pardiso control parameters. */
	int      iparm[64];
	double   dparm[64];
	int      solver;
	int      maxfct, mnum, phase, error, msglvl;

	/* Number of processors. */
	int      num_procs;

	/* Auxiliary variables. */
	char    *var;
	int      i, j, k;

	double   ddum;              /* Double dummy */
	int      idum;              /* Integer dummy. */
	
	error  = 0;
	solver = 0; /* use sparse direct solver */
	
	pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

	/* Numbers of processors, value of OMP_NUM_THREADS */
	var = getenv("OMP_NUM_THREADS");
	if(var != NULL)
		sscanf( var, "%d", &num_procs );
	else {
		printf("Set environment OMP_NUM_THREADS to 1\n");
		exit(1);
	}
	
	iparm[2]  = num_procs;

	iparm[10] = 1;      /* no scaling  */
	iparm[12] = 1;      /* no matching */
	maxfct = 1;         /* Maximum number of numerical factorizations.  */
	mnum   = 1;         /* Which factorization to use. */	
	msglvl = 0;         /* Print statistical information  */
	error  = 0;         /* Initialize error flag */

	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	int* perm = calloc (n,sizeof(int));
	for (i = 0; i < nrhs; ++i) 
		perm[i] = 1;
	iparm[30] = 1;
	
	double* aux_W = calloc(n*nrhs,sizeof(double));
	if (aux_W == NULL)	{ 
		printf ("Memory alloc failed. Exiting.. [PARDISO_top_tips]\n"); 
		exit(1);		
	}
	
	/* -------------------------------------------------------------------- */    
	/* ..  Reordering and Symbolic Factorization.  This step also allocates */
	/*     all memory that is necessary for the factorization.              */
	/* -------------------------------------------------------------------- */ 	
	phase = 11; 
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,perm,&nrhs,iparm,&msglvl,&ddum,&ddum,&error,dparm);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	
	/* -------------------------------------------------------------------- */    
	/* ..  Numerical factorization.                                         */
	/* -------------------------------------------------------------------- */    
	phase = 22;
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,perm,&nrhs,iparm,&msglvl,&ddum,&ddum,&error,dparm);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	
	/* -------------------------------------------------------------------- */    
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */ 		
	phase = 33;
	iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,perm,&nrhs,iparm,&msglvl,C,aux_W,&error,dparm);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	/* -------------------------------------------------------------------- */    
	/*  Arrange solution                                                    */
	/* -------------------------------------------------------------------- */  
	k = 0;
	for (i = 0 ; i < nrhs; i++)
		for (j = 0; j < nrhs; j++) 
			W[k++] = aux_W[i + j*n];
	
	free(aux_W);
	free(perm);

	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] -= 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}

	/* -------------------------------------------------------------------- */    
	/* ..  Termination and release of memory.                               */
	/* -------------------------------------------------------------------- */ 
	phase = -1;
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,&ddum,ia,ja,&idum,&nrhs,iparm,&msglvl,&ddum,&ddum,&error,dparm);
}

void PARDISO_numerical_factorization (MAT* A, void* pt[64], int iparm[64], double dparm[64], int nrhs)
{
	int      n = A->n;
	int    *ia = A->IA;
	int    *ja = A->JA;
	double  *a = A->AA;
	int    nnz = A->nz;

	int      i;
	double   ddum;               /* Double dummy */
	int      idum;               /* Integer dummy. */
	int      phase;
	int      mtype  = 11;        /* Real unsymmetric matrix */
	int      solver = 0;
	int      maxfct = 1;
	int      mnum   = 1;
	int      error  = 0;
	int      msglvl = 0;
	
	int         num_procs = atoi(getenv("OMP_NUM_THREADS"));
	iparm[2]  = num_procs;
	iparm[10] = 1;      /* no scaling  */
	iparm[12] = 1;      /* no matching */
	
	pardisoinit (pt, &mtype,&solver,iparm,dparm,&error);
	
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	/* -------------------------------------------------------------------- */    
	/* ..  Reordering and Symbolic Factorization.  This step also allocates */
	/*     all memory that is necessary for the factorization.              */
	/* -------------------------------------------------------------------- */ 	
	phase = 11; 
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,&idum,&nrhs,iparm,&msglvl,&ddum,&ddum,&error,dparm);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}		
	/* -------------------------------------------------------------------- */    
	/* ..  Numerical factorization.                                         */
	/* -------------------------------------------------------------------- */    
	phase = 22;
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,&idum,&nrhs,iparm,&msglvl,&ddum,&ddum,&error,dparm);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] -= 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
}

void PARDISO_back_substitution (MAT* A, void* pt[64], int iparm[64], double dparm[64], double* f, double* g, int nrhs)
{
	int      n = A->n;
	int    *ia = A->IA;
	int    *ja = A->JA;
	double  *a = A->AA;
	int    nnz = A->nz;

	int      i;
	int      idum;               /* Integer dummy. */
	int      phase;
	int      mtype  = 11;        /* Real unsymmetric matrix */
	int      maxfct = 1;
	int      mnum   = 1;
	int      error  = 0;
	int      msglvl = 0;
	
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	/* -------------------------------------------------------------------- */    
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */ 		
	phase = 33;
	iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,&idum,&nrhs,iparm,&msglvl,f,g,&error,dparm);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] -= 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
}


void PARDISO_back_substitution_top (MAT* A, void* pt[64], int iparm[64], double dparm[64], double* f, double* g, int nrhs, int k)
{
	int      n = A->n;
	int    *ia = A->IA;
	int    *ja = A->JA;
	double  *a = A->AA;
	int    nnz = A->nz;

	int      i;
	int      phase;
	int      mtype  = 11;        /* Real unsymmetric matrix */
	int      maxfct = 1;
	int      mnum   = 1;
	int      error  = 0;
	int      msglvl = 0;
	
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	int* perm = calloc (n,sizeof(int));
	for (i = 0; i < k; ++i) 
		perm[i] = 1;
	iparm[30] = 1;
	
	/* -------------------------------------------------------------------- */    
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */ 		
	phase = 33;
	iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,perm,&nrhs,iparm,&msglvl,f,g,&error,dparm);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] -= 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
}

void PARDISO_back_substitution_bottom (MAT* A, void* pt[64], int iparm[64], double dparm[64], double* f, double* g, int nrhs, int k)
{
	int      n = A->n;
	int    *ia = A->IA;
	int    *ja = A->JA;
	double  *a = A->AA;
	int    nnz = A->nz;

	int      i;
	int      phase;
	int      mtype  = 11;        /* Real unsymmetric matrix */
	int      maxfct = 1;
	int      mnum   = 1;
	int      error  = 0;
	int      msglvl = 0;
	
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	int* perm = calloc (n,sizeof(int));
	for (i = n - 1; i >= n - k; --i) 
		perm[i] = 1;
	iparm[30] = 1;
	
	/* -------------------------------------------------------------------- */    
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */ 		
	phase = 33;
	iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,a,ia,ja,perm,&nrhs,iparm,&msglvl,f,g,&error,dparm);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	for (i = 0; i < n; i++) 
		fprintf (stderr,"%lf\n",g[i]);
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */ 
	for (i = 0; i < n+1; i++) {
		ia[i] -= 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
}

void PARDISO_release_memory (MAT* A, void* pt[64], int iparm[64], double dparm[64], int nrhs)
{
	int      n = A->n;
	int    *ia = A->IA;
	int    *ja = A->JA;

	double   ddum;              /* Double dummy */
        int      idum;              /* Integer dummy. */

	int      phase;
	int      mtype  = 11;        /* Real unsymmetric matrix */
	int      maxfct = 1;
	int      mnum   = 1;
	int      error  = 0;
	int      msglvl = 0;
	
	/* -------------------------------------------------------------------- */    
	/* ..  Termination and release of memory.                               */
	/* -------------------------------------------------------------------- */ 
	phase = -1;
	pardiso (pt,&maxfct,&mnum,&mtype,&phase,&n,&ddum,ia,ja,&idum,&nrhs,iparm,&msglvl,&ddum,&ddum,&error,dparm);
}


