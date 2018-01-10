#include "protos.h"

/*----------------------------------------------------------------------------
 * Chains-on-chains partitioning
 *--------------------------------------------------------------------------*/
void PARTITIONING_MIN_MAX (int n, int *IA, int *JA, int nP, int** Fs, int** FL) 
{
	/* -------------------------------------------------------------------------------------
	 w: vetor de n pesos (nnz por linha)
	 s: vetor de k+1 separadores (cada particao i=0,k-1 possui as linhas s[i] a s[i+1]-1)
	 L: vetor de k cargas (peso total de cada particao)
	--------------------------------------------------------------------------------------*/
	
	int  i, j, B;
	int* w = (int *) malloc(n*sizeof(int));
	int* s = (int *) malloc((nP+1)*sizeof(int));
	int* L = (int *) malloc(nP*sizeof(int));
		
	for (i = 0; i < n; ++i)
		w[i] = IA[i+1] - IA[i];
	
	for (i = 0; i < nP; ++i) 
		s[i] = i;
	
	s[nP] = n;
	for (i = 0; i < nP; ++i)
	{
		L[i] = 0;
		for (j = s[i]; j < s[i+1]; ++j)
			L[i] += w[j];
	}
	do {
		for (i = 0, j = 1; j < nP; ++j) 
			if (L[j] > L[i]) i = j;
		B = L[i];
		j = i;
		
		if (s[i+1] == s[i]+1) break;
		
		while (L[i] >= B && i > 0)
		{
			L[i]   -= w[s[i]];
			L[i-1] += w[s[i]];
			s[i]    = s[i]+1;
			if (L[i] < B) i = i-1;
		}
	} while (L[0] < B);
	while (i != j)
	{
		if (L[i+1] < B) i = i+1;
		s[i]    = s[i]-1;
		L[i]   += w[s[i]];
		L[i-1] -= w[s[i]];
	}
	
	(*Fs) = s;
	(*FL) = L;
	free(w); 
	
	
	int outside = PARTITIONING_nnzout(IA, JA,s,nP);
		
	
	if (outside != 0)
	{
		fprintf(stderr, "\n Error. %d elements outside of blocks (must be zero). Exiting.. [PARTITIONING_MIN_MAX]\n\n",outside);
		exit(0);
	}
	
	
	return;
}

/*----------------------------------------------------------------------------
 * Verify if there are elements outside de diagonal or coupling blocks
 *--------------------------------------------------------------------------*/
int PARTITIONING_nnzout (int *IA, int *JA, int *s, int nP)
{
	int i, ii, j, c=0;

	if (nP == 1)
		return 0;

	for (i = s[0]; i < s[1]; ++i) 
	{ 
		for (j = IA[i]; j < IA[i+1]; ++j)
			if (JA[j]>=s[2]) c++; 
	}
	for (i = 1; i < nP-1; ++i)
	{
		for (ii = s[i]; ii < s[i+1]; ++ii)
		{
			for (j = IA[ii]; j < IA[ii+1]; ++j)
				if (JA[j] < s[i-1] || JA[j] >= s[i+2]) c++;
		}
	}
	for (i = s[nP-1]; i < s[nP]; ++i)
	{ 
		for (j = IA[i]; j < IA[i+1]; ++j) 
			if (JA[j]<s[nP-2]) c++;		
	}
	
	return c;
}



