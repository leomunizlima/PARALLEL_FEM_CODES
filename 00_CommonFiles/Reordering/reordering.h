#ifdef SSTranspEquation2D
	#include "../../01_SS_Transport_Equation/2D/01_CommonFiles/SSTranspEquation.h"
#endif
#ifdef TranspEquation2D
	#include "../../02_Transport_Equation/2D/01_CommonFiles/TranspEquation.h"
#endif
#ifdef EulerEquations2D
	#include "../../03_Euler_Equations/2D/01_CommonFiles/EulerEquations.h"
#endif
#ifdef SSNavierStokesEquations2D
	#include "../../04_SS_NavierStokes_Equations/2D/01_CommonFiles/SSNavierStokesEquations.h"
#endif

#include "symrcm.h"

/*---------------------------------------------------------------------------
 * ARRAY STRUCTURE
 *--------------------------------------------------------------------------*/
#ifndef ARRAY_H
#define ARRAY_H

typedef struct 
{
	int arr1;
	int arr2;
	int arr3;
} ARRAY;

typedef struct 
{
	double arr1;
	int    arr2;
} ARRAY2;

#endif /* ARRAY_H */

void mc73_fiedler(int* n, int* lirn, int* irn, int* ip, int* list, double* fvector, int* info, double* a);
void mc60ad_(int *n,int *cont1,int *irn, int *icptr,int *icntl,int *iw,int *info);
void mc60bd_(int *n, int *lirn, int *ja, int *ia, int *nsup, int *svar, int *vars, int *iw); 
void mc60cd_(int *n, int *nsup, int *lirn, int *ja, int *ia,int *list,int *jcntl, int *permsv, double *weight, int **pair, int *info, int *iw, double *w);
void mc60dd_(int *n, int *nsup, int *svar, int *vars,int *permsv, int *p, int *possv);
void mc60hd_(int *n, int *nsup, int *lirn, int *irn, int *icptr, int *vars, int *mask, int *ls, int *xls, int *list, int info[6]);
void reordering(ParametersType *Parameters, int *JA, int *IA, int *perm);
void REORDERING_SPECTRAL (int n, int nz, int *ja, int *ia, int *p);
void REORDERING_SYMRCM(int n, int nz,int *ja, int *ia, int *p);
int COMPARE_eig (const void * a, const void * b);
int COMPARE_array (const void * a, const void * b);
void MATRIX_ROW_permutation (int n, int nz, int *JA, int *IA, int *p, int *pT);
void MATRIX_SPARSE_ROW_permutation (int n, int nz, int *JA, int *IA, int *aux, int *pT);
void MATRIX_COL_permutation (int n, int nz, int *JA, int *IA, int *p, int *pT);
int MATRIX_bandwidth (int n, int *JA, int *IA);



