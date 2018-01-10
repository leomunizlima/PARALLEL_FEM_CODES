#ifndef SSTranspEquation_h
#define SSTranspEquation_h
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"

#define NDIM 2            // number related to the dimension of the problem (2: two-dimensional, 3: three-dimensional)
#define NNOEL 3           // number of nodes per element
#define NDOF 1            // number of degrees of freedom
#define PI 3.14159265359

typedef struct
{
	double x,y;	// x and y coordenations
	int Type;	//Type of the node: 0 means boundary, 1 means unknown, or >1 means the number of interfaces between MPI ranks
	int Send;	//Identify if a value needs to be send to other RANK (Send=1 No send, Send=2 send to RANK-1, Send=3 send to RANK+1)
	int invP_id; 	//Identification of the node: unknown, boundary, or interface between MPI ranks according to inverse permutation of the node in 
			// the soluction vector (only the unknowns have a valid number) 
}NodeType;

typedef struct
{
	int Vertex[NDIM+1];
	int Type;
}ElementType;

struct Node_List{
	int J; // vertice representando a cabeca do arco
	double value;
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
	char ProblemTitle[200];				//Problem title 
	char Solver[200];				//Solver 
	char Preconditioner[200];			//Preconditioner
	char Scaling[200];
	char reordering[200];				//Reordering scheme
	char MatrixVectorProductScheme[200];		//Matrix vector product: EBE, EDE, or CSR
	char StabilizationForm[200];			//Stabilization form: SUPG or DD
	char ShockCapture[200];				//ShockCapture: CAU or YZbeta
	char h_Shock[200];				//Many options of h parameters used in Shock Capture
	double SolverTolerance; 			//Tolerence of the solver
	double NonLinearTolerance;			//Tolerance of the non linear system
	int KrylovBasisVectorsQuantity;			//Number of vector used in Krylov Subspace
	int nnodes, nNodes;				//Number of nodes in the fem submesh of the MPI rank and the Global number of nodes
	int nel, nEl;					//Number of elements of the fem submesh of the MPI rank and the Global number of elements
	int neq, neq_bef, neqaux_bef, neqaux_aft;	//Number of equations of the linear system of the current MPI rank and previous MPI rank
	int nsend_bef, nsend_aft, nrecv_bef, nrecv_aft; //number of values to send and receive to RANK-1 and RANK+1
	int NEQ; 					//Number of equations in subsequent partitions (neq_bef + neq + neq_aft)
	int nedge;					//Number of edges of the fem submesh of the MPI rank
	int nnzero;					//Number of nonzeros coeffients in the local matrix associated with RANK
	int nnzero_bef;					//Number of nonzeros coeffients in the local matrix associated with RANK-1 
	int nnzero_aft;					//Number of nonzeros coeffients in the local matrix associated with RANK+1 
	int LinearMaxIter;				//Maximum of iterations of the linear system 
	int iterations;					//current number of the solver iterations
	int k_SPIKE_block_size;				//Size of coupling blocks of SPIKE preconditioner
	int bandwidth_bef, bandwidth_aft;		//Half bandwidth before and after reordering
	int RANK; 					// ID of MPI rank
	int NPROC; 					// Number of MPI ranks
	char **arguments;
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
	int **Scheme_by_Element;	//Used in EDE or CSR format
	int **Scheme_by_Element_bef;	//Used in EDE or CSR format (associated with RANK-1)
	int **Scheme_by_Element_aft;	//Used in EDE or CSR format (associated with RANK+1)
	int **order;			//Used in EDE format
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
	int *Perm;			//Used in CSR format to ensure matrix permutation or local reordering to reduce cache miss
	int *invPerm;			//Used in local reordering to reduce cache miss
	double **fullC;			//Block before of SPIKE Preconditioner (size k x k)
	double *fullCaux;			//Aux allocation to C
	double **fullB;			//Block after of SPIKE Preconditioner (size k x k)
	double *fullBaux;			//Aux allocation to B
	double **LUe;
	double *LUeaux;
	MAT *A_SPIKE, *B_SPIKE, 
	    *C_SPIKE, *S_SPIKE;		//SPIKE structures
}MatrixDataType;

typedef struct
{
	int **lm;
	int *lmaux;
	NodeType *Node;
	ElementType *Element;
	double *delta_old; //store value of norm of resudium in previous iterarion for each element 
	double *Re_old;
	double *F; // right hand side vector with size of the solution plus RANK before and after
	double *U; // vector of the solution with size of the solution plus RANK before and after
	double *u; // vector of the solution (current RANK)
	double *f; // right hand side vector (current RANK)
	double *SendBuffer_bef, *SendBuffer_aft, *RecvBuffer_bef, *RecvBuffer_aft; //Buffers to send and receive data to RANK-1 and RANK+1
	int *IdSend_bef, *IdSend_aft, *IdRecv_bef, *IdRecv_aft;//Vector to identify positions Node[I].invP_id in the buffers
}FemStructsType;

typedef struct
{
	double (*Condutivity)(void);
	double (*Font)(double, double, double, double, double, double); 
	double (*Reaction)(void);
	void (*Velocity)(double, double, double []);
	double (*upresc)(double, double);
	double (*h_shock_calculation)(double, double, double, double, double, double, double, double, double, double, double, double);
	double (*ShockCapture)(double *, double *, int, double,  double, double, double, double, double, double, double, double, double, double, 
			       double, double, double, double, double, double, double);
	void (*assembly)(ParametersType *, MatrixDataType *, FemStructsType *, int, double [3][3]);
	int (*mv)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precondR)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond_setup)(ParametersType *Parameters, MatrixDataType *, FemStructsType *, int, double *);
	int (*scaling)(ParametersType *, MatrixDataType *, FemStructsType *);
	int (*unscaling)(ParametersType *, MatrixDataType *, FemStructsType *, double *);
}FemFunctionsType;

typedef struct
{
	int (*solver) (ParametersType *, MatrixDataType *, FemStructsType*, FemFunctionsType*, double *, double *); 
	int (*Build_K_F)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);
}FemOtherFunctionsType;


int Preprocess(int, int, int, char **, ParametersType **, MatrixDataType **, FemStructsType **, FemFunctionsType **, FemOtherFunctionsType **);

int Process(int, int, ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Postprocess(int, int, ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Paraview_Output(int, int, ParametersType *, FemStructsType *, FemFunctionsType *);

int Fill_LM(ParametersType *, int **, NodeType *, ElementType *);

void csr_Initialization(ParametersType *, MatrixDataType *, FemStructsType *);

void csr_List_insertA(NodeListType **, int, int, int *);

int csr_search(int, int, NodeListType *);

int setProblem(ParametersType *, FemFunctionsType *);

int setMatrixVectorProductType(ParametersType *, FemFunctionsType *);

int setSolver(ParametersType *, FemOtherFunctionsType *); 

int setPreconditioner(ParametersType *, FemFunctionsType *);

int setScaling(ParametersType *, FemFunctionsType *);

int setStabilizationForm(ParametersType *, FemFunctionsType *, FemOtherFunctionsType *); 

void ebe_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double [3][3]);

void csr_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double [3][3]);

void F_assembly(NodeType *, int, int, int, double, double, double, double *, double [3][3],double [3], double [3], FemFunctionsType *);

double CAU_ShockCapture(double *, double *, int, double, double, double, double, double, double, double, double, double, double, double, double, double, 
			double, double, double, double, double);

double CAU_DD_ShockCapture(double *, double *, int, double, double, double, double, double, double, double, double, double, double, double, double, double, 
			double, double, double, double, double);

double YZBeta_ShockCapture(double *, double *, int, double, double, double, double, double, double, double, double, double, double, double, double, double, 
			double, double, double, double, double);

double h_shock_2sqrtArea(double, double, double, double, double, double, double, double, double, double, double, double);

double h_shock_Option1(double, double, double, double, double, double, double, double, double, double, double, double);

double h_shock_Option2(double, double, double, double, double, double, double, double, double, double, double, double);

int setzeros(ParametersType *, MatrixDataType *);

int Build_K_F_SUPG(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_K_F_DD(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

void MPI_VectorUpdate(int, int, ParametersType *, FemStructsType *, double *, int);

int COMPARE_array_in_Preprocess (const void * a, const void * b);

void Permutation_of_LM(ParametersType *, FemStructsType *, int *, int *);

#endif




