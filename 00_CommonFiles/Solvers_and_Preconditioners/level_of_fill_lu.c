#include "preconditioners.h"
#include "ilup.h"
#include "../Allocation_Operations/allocations.h"

int level_of_fill_lu (ParametersType *Parameters, MatrixDataType *MatrixData)
{
	int n = Parameters->neq;
	int* levls = NULL;
	int* jbuf  = NULL;
	int* iw    = NULL;
	int** ulvl;  /*  stores lev-fils for U part of ILU factorization*/
	SparMAT *L, *U, *A; 
	int i, j, k, col, ip, it, jpiv;
	int incl, incu, jmin, kmin; 

  
	levls = (int*)  mycalloc("levls in level_of_fill_lu", n, sizeof(int));
	jbuf  = (int*)  mycalloc("jbuf in level_of_fill_lu", n,sizeof(int)); 
	ulvl  = (int**) mycalloc("ulvl in level_of_fill_lu", n,sizeof(int*));
	iw = (int*)  mycalloc(n*sizeof(int));
	

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


	MatrixData->ILUp->work = iw;

	return 0;
}  

