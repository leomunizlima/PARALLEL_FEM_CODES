#ifndef SSNavierStokesEquations_h
#define SSNavierStokesEquations_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include <time.h>
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"

#define NDIM 2            // number related to the dimension of the problem (2: two-dimensional, 3: three-dimensional)
#define NNOEL 3           // number of nodes per element
#define NDOF 3            // number of degrees of freedom
#define PI 3.14159265359

typedef struct
{
	double x,y;    // coordinates of the node
	int Type[NDOF];           // mark of the node according to boundary conditions
	int Send[NDOF];	//Identify if a value needs to be send to other RANK (Send=1 No send, Send=2 send to RANK-1, Send=3 send to RANK+1)
	int invP_id[NDOF];  // vector that identifies whether the node is prescribed or not for each property
                         // 0: velocity in the x direction
                         // 1: velocity in the y direction
                         // 2: pression
}NodeType;             

typedef struct
{
	int Vertex[NNOEL];  // holds the global number of nodes that make up the element
	int Type;           // mark of the element
}ElementType;

struct Node_List{
	double value;
	int J; // vertice representando a cabeca do arco
	struct Node_List *next;
};

typedef struct Node_List NodeListType;

typedef struct{
	int array1;
	int array2;
}ARRAYType;

typedef struct{
	/* Pardiso parameters. */
	void*     A_pt[64];
	int       A_iparm[64];
	double    A_dparm[64];
	void*     S_pt[64];
	int       S_iparm[64];
	double	  S_dparm[64];
}PardisoVariablesType;

typedef struct
{
	char ProblemTitle[200];                // Problem name to be solved
	char Solver[200];                      // type of method used in the solution
	char MatrixVectorProductScheme[200];   // the global matrix storage form
	char StabilizationForm[200];           // type of stabilization method
	char Preconditioner[200];           // preconditioners: yes - use or not - don't use
	char Experiments[200];
	char reordering[200];			// Reordering for CSR (NOT, Spectral, Weigthed Spectral (WSO) or RCM)
	double SolverTolerance;                // tolerance for the solution method
	double NonLinearTolerance;            // tolerance for the loop of correction
	double ReynoldsNumber; 			//Reynold's number
	int KrylovBasisVectorsQuantity;        // Krylov number of vectors in the basis for the restart
	int nnodes, nNodes;                            // nnodes: number of nodes/ nNodes: global number of nodes 
	int nel, nEl;                               // nel: number of element / nEl: global number of elements
	int neq, neq_bef, neqaux_bef, neqaux_aft;	//Number of equations of the linear system of the current MPI rank and previous MPI rank
	int nsend_bef, nsend_aft, nrecv_bef, nrecv_aft; //number of values to send and receive to RANK-1 and RANK+1
	int NEQ; 					//Number of equations in subsequent partitions (neq_bef + neq + neq_aft)
	int nnzero;					//Number of nonzeros coeffients in the local matrix associated with RANK
	int nnzero_bef;					//Number of nonzeros coeffients in the local matrix associated with RANK-1 
	int nnzero_aft;					//Number of nonzeros coeffients in the local matrix associated with RANK+1 
	int iterations;                        // iterations: total number of iteration 
	int k_SPIKE_block_size;				//Size of coupling blocks of SPIKE preconditioner
	int RANK; 					// ID of MPI rank
	int NPROC; 					// Number of MPI ranks
	int LinearMaxIter;                           // itermax: maximum number of iteration
	int NonLinearMaxIter;                           // itermax: maximum number of iteration
	int bandwidth_bef, bandwidth_aft;	// half bandwidth before and after reordering
	PardisoVariablesType *PardisoVariables;		//Variables to be used inside Pardiso solver
}ParametersType;

typedef struct
{
	double **A;   			//Used in EBE or EDE format
	double *Aaux;			//Used in EBE or EDE format. It ensures allocations in line 
	double *AA;			//Used only in CSR format
	double *AA_bef;			//Used only in CSR format (associated with RANK-1)
	double *AA_aft;			//Used only in CSR format (associated with RANK+1)
	double *Diag; 			//Used in diagonal preconditioning
	double *invDiag;		//Used in diagonal preconditioning
	double **BlockDiag;
	double *BlockDiagAux;
	double **invBlockDiag;
	double *invBlockDiagAux;
	double **invDe;
	double *invDeaux;
	double **LUe;
	double *LUeaux;
	int **Id;
	int *IdAux;
	int **Scheme_by_Element;	//Used in EDE or CSR format
	int **Scheme_by_Element_bef;	//Used in EDE or CSR format (associated with RANK-1)
	int **Scheme_by_Element_aft;	//Used in EDE or CSR format (associated with RANK+1)
	int **order;
	int *JA;			//Used only in CSR format
	int *IA;			//Used only in CSR format
	int *JA_bef;			//Used only in CSR format (associated with RANK-1)
	int *IA_bef;			//Used only in CSR format (associated with RANK-1)
	int *JA_aft;			//Used only in CSR format (associated with RANK+1)
	int *IA_aft;			//Used only in CSR format (associated with RANK+1)
	int *aux_bef;			//Used only in CSR format (associated with RANK-1)
	int *aux_aft;			//Used only in CSR format (associated with RANK+1)
	SparILU *ILUp;			//Used only in CSR format to ILUp preconditioning 
	SparMAT *mat;
	MAT *Ailu;
	int *Perm;
	int *invPerm;
	double **fullC;			//Block before of SPIKE Preconditioner (size k x k)
	double *fullCaux;			//Aux allocation to C
	double **fullB;			//Block after of SPIKE Preconditioner (size k x k)
	double *fullBaux;			//Aux allocation to B
	MAT *A_SPIKE, *B_SPIKE, 
	    *C_SPIKE, *S_SPIKE;		//SPIKE structures
}MatrixDataType;

typedef struct{
	double **M2, **R2, *invN2, *delta_old_NMV;	
	double tolerance;
}AuxBuildStructuresType;

typedef struct
{
	int **lm;
	int *lmaux;
	NodeType *Node;
	ElementType *Element;
	double *F; // right hand side vector with size of the solution plus RANK before and after
	double *U; // vector of the solution with size of the solution plus RANK before and after
	double *u; // vector of the solution (current RANK)
	double *f; // right hand side vector (current RANK)
	double *uB;
	double *delta_old;
	AuxBuildStructuresType *AuxBuild;
	double *SendBuffer_bef[2], *SendBuffer_aft[2], *RecvBuffer_bef[2], *RecvBuffer_aft[2]; //Buffers to send and receive data to RANK-1 and RANK+1
	int *IdSend_bef, *IdSend_aft, *IdRecv_bef, *IdRecv_aft;//Vector to identify positions Node[I].invP_id in the buffers
}FemStructsType;

typedef struct
{
	double (*v1presc)(double, double);
	double (*v2presc)(double, double);
	double (*ppresc)(double, double);
	void (*assembly)(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[9]);
	int (*mv)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precondR)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond_setup)(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
}FemFunctionsType;

typedef struct
{
	int (*solver) (ParametersType *, MatrixDataType *, FemStructsType*, FemFunctionsType*, double *, double *);
	int (*Build)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);
}FemOtherFunctionsType;


int Preprocess(int, int, int, char **, ParametersType **, MatrixDataType **, FemStructsType **, FemFunctionsType **, FemOtherFunctionsType **);

int Process(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Postprocess(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Fill_LM(ParametersType *, int **, NodeType *, ElementType *);

int setProblem(ParametersType *, FemFunctionsType *);

int setMatrixVectorProductType(ParametersType *, FemFunctionsType *);

int setSolver(ParametersType *, FemOtherFunctionsType *);

int setPreconditioner(ParametersType *, FemFunctionsType *);

int setStabilizationForm(ParametersType *,FemFunctionsType *, FemOtherFunctionsType *);

int setzeros(ParametersType *, MatrixDataType *);

void csr_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[9]);

void csr_Initialization(ParametersType *, MatrixDataType *, FemStructsType *);

void csr_List_insertA(NodeListType **, int , int , int *);

int csr_search(int, int, NodeListType *);

void ebe_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[9]);

int Build_K_F_SUPG_PSPG(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_K_F_VMS_DCDD(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

void eval_U(ParametersType *,FemStructsType *, FemFunctionsType *, double *);

int Paraview_Output(int, int, ParametersType *, FemStructsType *, FemFunctionsType *);

int Paraview_Output_3D(int, int, ParametersType *, FemStructsType *, FemFunctionsType *);

int COMPARE_array_in_Preprocess (const void * a, const void * b);

void MPI_VectorUpdate(int, int, ParametersType *, FemStructsType *, double *, int);

#endif

