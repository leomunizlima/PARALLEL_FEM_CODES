#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "protos.h"

/*----------------------------------------------------------------------------
 * PARDISO
 *--------------------------------------------------------------------------*/
void     pardisoinit                     (void*,int*,int*,int*,double*,int*);
void     pardiso                         (void*,int*,int*,int*,int*,int*,double*,int*,int*,int*,int*,int*,int*,double*,double*,int*,double*);
void     pardiso_chkmatrix               (int*,int*,double*,int*,int*,int*);
void     pardiso_chkvec                  (int*,int*,double*,int*);
void     pardiso_printstats              (int*,int*,double*,int*,int*,int*,double*,int*);
void     read_matrix                     (double**,int**,int**,double**,int*,int*,int*,FILE*);
void     PARDISO_bottom_tips             (MAT* A, double* B, double* V, int nrhs);
void     PARDISO_top_tips                (MAT* A, double* C, double* W, int nrhs);
void     PARDISO_numerical_factorization (MAT* A, void* pt[64], int iparm[64], double dparm[64], int nrhs);
void     PARDISO_back_substitution       (MAT* A, void* pt[64], int iparm[64], double dparm[64], double* f, double* g, int nrhs);
void     PARDISO_release_memory          (MAT* A, void* pt[64], int iparm[64], double dparm[64], int nrhs);

