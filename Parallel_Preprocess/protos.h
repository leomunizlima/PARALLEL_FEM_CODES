#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "symrcm.h"

#define NDOF_MAX 4

typedef struct
{
	double x,y;
	int Type[NDOF_MAX];
	int id[NDOF_MAX];
}NodeType;

typedef struct
{
	int Vertex[3];
	int Type;
}ElementType;

typedef struct List
{
	int J;
	double value;
	struct List *next;
} ListType;

typedef struct 
{
	double arr1;
	int    arr2;
	int    arr3;
} ARRAY;

void REORDERING_SPECTRAL (int n, int nnz, int *ia, int *ja, int *P);
void REORDERING_SYMRCM(int n, int nnz, int *ja, int *ia, int *p);
void insertNode(ListType **L, int i, int j, int *nnz);
void mount_IA_JA(int n, int nnz, ListType **L,  int *IA, int *JA);
void MATRIX_ROW_permutation (int n, int nnz, int *IA, int *JA, int* p);
void MATRIX_COL_permutation (int n, int nnz, int *IA, int *JA, int* p);
int MATRIX_bandwidth (int n, int *IA, int *JA);
int COMPARE_array (const void * a, const void * b);
void PARTITIONING_MIN_MAX (int n, int *IA, int *JA, int nP, int** Fs, int** FL);
int PARTITIONING_nnzout (int *IA, int *JA, int *s, int nP);
int COUNTER(ListType **L, int n, int K);
int catch_Type(NodeType *Node, int I, int K, int JJ, int *d, int *invP); 
int find_neq_mesh(int *vec, int n);
int Search(ListType **L, int I, int K);
void free_List(ListType **L, int n);

