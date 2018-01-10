#ifndef PROTOS_H
#define PROTOS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define min(a,b) (((a)>(b))?(b):(a))
#define max(a,b) (((a)>(b))?(a):(b))

typedef struct
{
	int n, m, nz;
	double *D;
	double *AA;
	int *JA;
	int *IA;	
} MAT;

typedef struct
{
	int           n;
	int*    nzcount;  /* length of each row                          */
	int**        ja;  /* pointer-to-pointer to store column indices  */
	double**     ma;  /* pointer-to-pointer to store nonzero entries */
} SparMAT;

typedef struct
{
	int         n;
	SparMAT*    L;   /* L part elements   */
	double*     D;   /* diagonal elements */
	SparMAT*    U;   /* U part elements   */
	int*     work;   /* working buffer    */
} SparILU;

#endif /* PROTOS_H */

