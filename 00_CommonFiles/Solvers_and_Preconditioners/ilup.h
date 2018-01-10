/*----------------------------------------------------------------------------
 * MATRIX HEADER FILE
 *--------------------------------------------------------------------------*/
#ifndef PRECONDITIONERS_H
#define PRECONDITIONERS_H

#include "protos.h"

void SPARMAT_setup  (SparMAT* mat, int n);
void SPARILU_setup  (SparILU* lu, int n);
void SPARIUL_setup  (SparILU* ul, int n);
void SPARILU_row    (SparILU* lu, int nrow);
void SPARIUL_row    (SparILU* ul, int nrow);
void CSRto_SPARMAT (MAT* A, SparMAT* mat);
void SPARILU_toCSR (SparILU* lu, MAT* L, MAT* U);
int  LEVEL_OF_FILL_lu  (SparMAT* csmat, SparILU* lu, int p);
int  LEVEL_OF_FILL_ul  (SparMAT* csmat, SparILU* ul, int p);
void ILUP           (SparMAT* csmat, SparILU* lu, int p);
void IULP           (SparMAT* csmat, SparILU* lu, int p);
void SPARMAT_clean  (SparMAT* mat);
void SPARILU_clean  (SparILU* lu);
void SPARILU_print  (SparILU* lu);
void QSPLIT         (double *a, int *ind, int n, int Ncut);
void ILUT           (SparMAT* csmat, SparILU* lu, int lfil, double tol);

#endif /* PRECONDITIONERS_H */
