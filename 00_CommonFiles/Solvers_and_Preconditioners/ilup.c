/*----------------------------------------------------------------------------
 * ILUP PRECONDITIONER
 *--------------------------------------------------------------------------*/
#include "ilup.h"

/*----------------------------------------------------------------------------
 * Initialize SparMAT structs
 *--------------------------------------------------------------------------*/
void SPARMAT_setup (SparMAT* mat, int n)
{
	mat->n       = n;
	if ((mat->nzcount = (int *) calloc(n,sizeof(int)))==NULL) { printf("ERRO ALOCACAO of mat->nzcount!\n"); exit(1);}
	if ((mat->ja      = (int **) calloc(n,sizeof(int*)))==NULL) {printf ("ERRO ALOCACAO of mat->ja\n"); exit(1);}
	if ((mat->ma      = (double **) calloc(n,sizeof(double*)))==NULL) {printf ("ERRO ALOCACAO of mat->ma\n"); exit(1);}

}

/*----------------------------------------------------------------------------
 * Initialize SparILU structs
 *--------------------------------------------------------------------------*/
void SPARILU_setup (SparILU* lu, int n)
{
	lu->n    = n;
	lu->D    = calloc(n,sizeof(double));
	lu->L    = (SparMAT*) malloc(sizeof(SparMAT));
	SPARMAT_setup(lu->L, n);
	lu->U    = (SparMAT*) malloc(sizeof(SparMAT));
	SPARMAT_setup(lu->U, n);
	lu->work = calloc(n,sizeof(int));
}

/*----------------------------------------------------------------------------
 * Prepare space of a row according to the result of level structure
 *--------------------------------------------------------------------------*/
void SPARILU_row (SparILU* lu, int nrow)
{
	int nzcount;
	nzcount = lu->L->nzcount[nrow];
	if (nzcount){
		free(lu->L->ma[nrow]);
		lu->L->ma[nrow] = calloc(nzcount,sizeof(double));
	}
	nzcount = lu->U->nzcount[nrow];
	if (nzcount){
		free(lu->U->ma[nrow]);
		lu->U->ma[nrow] = calloc(nzcount,sizeof(double));
	}
}

/*----------------------------------------------------------------------------
 * Convert CSR matrix to SparMAT struct
 *--------------------------------------------------------------------------*/
void CSRto_SPARMAT_setup (MAT* A, SparMAT* mat)
{
	int i, j, j1, len;
	int n = A->n;
	double* bra;
	int*    bja;
	
	/* setup data structure for mat (SparMAT*) struct */
	SPARMAT_setup(mat, n);
  
	for (j = 0; j < n; ++j)
	{
		len = A->IA[j+1] - A->IA[j];
		mat->nzcount[j] = len;
		if (len > 0) 
		{
			bja = malloc(len*sizeof(int));
			bra = malloc(len*sizeof(double));
			i   = 0;
			for (j1 = A->IA[j]; j1 < A->IA[j+1]; ++j1)
			{
				bja[i] = A->JA[j1];
				bra[i] = A->AA[j1];
				i++;
			}
			mat->ja[j] = bja;
			mat->ma[j] = bra;
		}
	}    
}

/*----------------------------------------------------------------------------
 * Convert CSR matrix to SparMAT struct
 *--------------------------------------------------------------------------*/
void CSRto_SPARMAT (MAT* A, SparMAT* mat)
{
	int i, j, j1;
	int n = A->n;
	
	for (j = 0; j < n; ++j)
	{
		i   = 0;
		for (j1 = A->IA[j]; j1 < A->IA[j+1]; ++j1)
		{
			mat->ma[j][i] = A->AA[j1];
			i++;
		}
	}    
}


/*----------------------------------------------------------------------------
 * Convert SparILU struct to CSR
 *--------------------------------------------------------------------------*/
void SPARILU_toCSR  (SparILU* lu, MAT* L, MAT* U)
{
	int i, j, Lcount = 0, Ucount = 0;
		
	for (i = 0; i < lu->n; ++i)
	{
		for (j = 0; j < lu->L->nzcount[i]; ++j) Lcount++;
		for (j = 0; j < lu->U->nzcount[i]; ++j) Ucount++;
	}
		
	L->AA  = (double *) calloc(Lcount+lu->n, sizeof (double));
	L->D   = (double *) calloc(lu->n,        sizeof (double));
	L->JA  = (int    *) calloc(Lcount+lu->n, sizeof (int));
	L->IA  = (int    *) calloc(lu->n+1,      sizeof (int));

	U->AA  = (double *) calloc(Ucount+lu->n, sizeof (double));
	U->D   = (double *) calloc(lu->n,        sizeof (double));
	U->JA  = (int    *) calloc(Ucount+lu->n, sizeof (int));
	U->IA  = (int    *) calloc(lu->n+1,      sizeof (int));
	
	L->n     = lu->n;
	L->nz    = Lcount;
	Lcount   = 0;
	L->IA[0] = 0;
	for (i = 0; i < lu->n; ++i)
	{
		for (j = 0; j < lu->L->nzcount[i]; ++j)
		{
			L->AA[Lcount] = lu->L->ma[i][j];
			L->JA[Lcount] = lu->L->ja[i][j];
			Lcount++;
		}
		L->D[i]    = 1.0;
		L->IA[i+1] = Lcount;
	}
	
	U->n     = lu->n;
	U->nz    = Ucount;
	Ucount   = 0;
	U->IA[0] = 0;
	for (i = 0; i < lu->n; ++i)
	{
		for (j = 0; j < lu->U->nzcount[i]; ++j)
		{
			U->AA[Ucount] = lu->U->ma[i][j];
			U->JA[Ucount] = lu->U->ja[i][j];
			Ucount++;
		}
		U->D[i]    = 1/lu->D[i];
		U->IA[i+1] = Ucount;
	}	
}

/*----------------------------------------------------------------------------
 * symbolic ilu factorization to calculate structure of ilu matrix
 * for specified level of fill
	 *--------------------------------------------------------------------------*/
int LEVEL_OF_FILL_lu (SparMAT* csmat, SparILU* lu, int p)
{
	int n = csmat->n;
	int* levls = NULL;
	int* jbuf  = NULL;
	int* iw    = lu->work;
	int** ulvl;  /*  stores lev-fils for U part of ILU factorization*/
	SparMAT* L = lu->L; 
	SparMAT* U = lu->U;
	/*--------------------------------------------------------------------
	* n        = number of rows or columns in matrix
	* inc      = integer, count of nonzero(fillin) element of each row
	*            after symbolic factorization
	* ju       = entry of U part of each row
	* lvl      = buffer to store levels of each row
	* jbuf     = buffer to store column index of each row
	* iw       = work array
	*------------------------------------------------------------------*/
	int i, j, k, col, ip, it, jpiv;
	int incl, incu, jmin, kmin; 
  
	levls = (int*)  malloc(n*sizeof(int));
	jbuf  = (int*)  malloc(n*sizeof(int)); 
	ulvl  = (int**) malloc(n*sizeof(int*));

	/* initilize iw */
	for(j = 0; j < n; j++) iw[j] = -1;
	for(i = 0; i < n; i++) 
	{
		incl = 0;
		incu = i; 
		/*-------------------- assign lof = 0 for matrix elements */
		for(j = 0; j < csmat->nzcount[i]; j++) 
		{
			col = csmat->ja[i][j];
			if(col < i) 
			{
				/*-------------------- L-part  */
				jbuf[incl]  = col;
				levls[incl] = 0;
				iw[col]     = incl++;
			} 
			else if (col > i) 
			{ 
				/*-------------------- U-part  */
				jbuf[incu]  = col;
				levls[incu] = 0;
				iw[col]     = incu++;
			} 
		}
		/*-------------------- symbolic k,i,j Gaussian elimination  */ 
		jpiv = -1; 
		while (++jpiv < incl)
		{
			k = jbuf[jpiv] ; 
			/*-------------------- select leftmost pivot */
			kmin = k;
			jmin = jpiv; 
			for(j = jpiv + 1; j< incl; j++)
			{
				if(jbuf[j] < kmin)
				{
					kmin = jbuf[j];
					jmin = j;
				}
			}
			/*-------------------- swap  */  
			if(jmin != jpiv)
			{
				jbuf[jpiv]  = kmin; 
				jbuf[jmin]  = k; 
				iw[kmin]    = jpiv;
				iw[k]       = jmin; 
				j           = levls[jpiv] ;
				levls[jpiv] = levls[jmin];
				levls[jmin] = j;
				k           = kmin; 
			}
			/*-------------------- symbolic linear combinaiton of rows  */
			for(j = 0; j < U->nzcount[k]; j++)
			{
				col = U->ja[k][j];
				it  = ulvl[k][j]+levls[jpiv]+1 ; 
				if(it > p) continue; 
				ip  = iw[col];
				if( ip == -1 )
				{
					if(col < i)
					{
						jbuf[incl]  = col;
						levls[incl] = it;
						iw[col]     = incl++;
					} 
					else if(col > i)
					{
						jbuf[incu]  = col;
						levls[incu] = it;
						iw[col]     = incu++;
					} 
				}
				else
					levls[ip] = min(levls[ip], it); 
			}
		}   /* end - while loop */
		/*-------------------- reset iw */
		for(j = 0; j < incl; j++) iw[jbuf[j]] = -1;
		for(j = i; j < incu; j++) iw[jbuf[j]] = -1;
		/*-------------------- copy L-part */ 
		L->nzcount[i] = incl;
		if(incl > 0 )
		{
			L->ja[i] = (int *) malloc(incl*sizeof(int));
			memcpy(L->ja[i], jbuf, incl*sizeof(int));
		}
		/*-------------------- copy U - part        */ 
		k = incu-i; 
		U->nzcount[i] = k; 
		if( k > 0 )
		{
			U->ja[i] = (int *) malloc(k*sizeof(int));
			memcpy(U->ja[i], jbuf+i, k*sizeof(int));
			/*-------------------- update matrix of levels */
			ulvl[i]  = (int *) malloc(k*sizeof(int)); 
			memcpy(ulvl[i], levls+i, k*sizeof(int));
		}
	}
  
	/*-------------------- free temp space and leave --*/
	free(levls);
	free(jbuf);
	for(i = 0; i < n-1; i++)
	{
		if (U->nzcount[i])
			free(ulvl[i]) ; 
	}
	free(ulvl); 

	return 0;
}  


/*----------------------------------------------------------------------------
 * Free up memory allocated for SpaFmt structs.
 *--------------------------------------------------------------------------*/
void SPARMAT_clean (SparMAT* mat)
{
	int i;
	if (mat == NULL)
	{ 
		free(mat); 
		return; 
	}
	if (mat->n < 1)
	{ 
		free(mat);
		return; 
	}
	for (i = 0; i < mat->n; ++i) {
		if (mat->nzcount[i] > 0)
		{
			if (mat->ma) free(mat->ma[i]);
			free(mat->ja[i]);
		}
	}    
	if (mat->ma) free(mat->ma);
	free(mat->ja);
	free(mat->nzcount);
	free(mat);
	return;
}

/*----------------------------------------------------------------------------
 * Free up memory allocated for SparILU structs.
 *--------------------------------------------------------------------------*/
void SPARILU_clean (SparILU* lu)
{
	if (lu == NULL)
	{ 
		free(lu);
		return; 
	}
	if (lu->n < 1)
	{ 
		free(lu); 
		return; 
	}
	if (lu->D)    free(lu->D);
	SPARMAT_clean (lu->L);
	SPARMAT_clean (lu->U);
	if (lu->work) free(lu->work);
	free (lu);
	return;
}


	
