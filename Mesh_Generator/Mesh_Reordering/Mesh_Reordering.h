#include<stdlib.h>
#include<stdio.h>
#include<string.h>

typedef struct
{
	int Vertex[3];
	int Type;
}ElementType;

FILE *myfopen(char *, char *);
void *mycalloc(char *, int, int);
void reordering_rcm(int, int, int, char *, double *, double *, int **, ElementType *);
void reading_mesh(char *, int, int, int, double *, double *, int **, ElementType *);
void writing_mesh(int, int, int, int **, ElementType *, double *, double *, char *);
void checking_mesh(int *, int, double *, double *, int **, ElementType *);


