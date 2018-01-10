#ifndef EulerEquations_h
#define EulerEquations_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include <time.h>
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"

#define NDIM 2            // number related to the dimension of the problem (2: two-dimensional, 3: three-dimensional)
#define NNOEL 3           // number of nodes per element
#define NDOF 4            // number of degrees of freedom
#define PI 3.14159265359

typedef struct
{
	double x,y;    // coordinates of the node
	int Type[NDOF];           // mark of the node according to boundary conditions
	int Send[NDOF];	//Identify if a value needs to be send to other RANK (Send=1 No send, Send=2 send to RANK-1, Send=3 send to RANK+1)
	int invP_id[NDOF];  // vector that identifies whether the node is prescribed or not for each property
                         // 0: density
                         // 1: velocity in the x direction
                         // 2: velocity in the y direction
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
	char TimeIntegration[200];             // Time integration method
	char MatrixVectorProductScheme[200];   // the global matrix storage form
	char StabilizationForm[200];           // type of stabilization method
	char ShockCapture[200];            // type discontinuities capture Operator
	char Preconditioner[200];           // preconditioners: yes - use or not - don't use
	char Scaling[200];
	char Experiments[200];
	char StopMulticorrection[200];         // Fala se o loop espacial para pela norma ou por um numero fixo de iteracao - NORM: para pela norma; ITERATION: para pela iteracao
	char reordering[200];			// Reordering for CSR (NOT, Spectral, Weigthed Spectral (WSO) or RCM)
	char Dimensionless[200];		// Determines whether or not a problem is dimensionless
	char StopAtSteadyState[200];		// YES or NO if you want to stop in steady state or final time
	double SolverTolerance;                // tolerance for the solution method
	double NonLinearTolerance;            // tolerance for the loop of correction
	double TimeIntegrationTolerance;       // tolerance for the loop of time integration
	double StabilizationTolerance;           // tolerance for coefficient used in the stabilization form
	double Mach;				//Mach number of reference
	double Alpha;                          // parameter determining the stability control and accuracy time integration method
	double Alpha_Build;					// parameter Alpha to be used inside Build functions
	double DeltaT;                         // time step
	double DeltaT_Build;					// time step to be used inside Build functions
	double FinalTime;                         // final time
	double CurrentTime;                         // current time (to be used in steady state situation)
	double invY[4];                        // vector that stores 1/U used in the capture operator YZBeta
	int KrylovBasisVectorsQuantity;        // Krylov number of vectors in the basis for the restart
	int nnodes, nNodes;                            // nnodes: number of nodes/ nNodes: global number of nodes 
	int nel, nEl;                               // nel: number of element / nEl: global number of elements
	int neq, neq_bef, neqaux_bef, neqaux_aft;	//Number of equations of the linear system of the current MPI rank and previous MPI rank
	int nsend_bef, nsend_aft, nrecv_bef, nrecv_aft; //number of values to send and receive to RANK-1 and RANK+1
	int NEQ; 					//Number of equations in subsequent partitions (neq_bef + neq + neq_aft)
	int neqrho;                            // neqrho: number of equations relating the density
	int nnzero;					//Number of nonzeros coeffients in the local matrix associated with RANK
	int nnzero_bef;					//Number of nonzeros coeffients in the local matrix associated with RANK-1 
	int nnzero_aft;					//Number of nonzeros coeffients in the local matrix associated with RANK+1 
	int iterations;                        // iterations: total number of iteration 
	int k_SPIKE_block_size;				//Size of coupling blocks of SPIKE preconditioner
	int RANK; 					// ID of MPI rank
	int NPROC; 					// Number of MPI ranks
	int LinearMaxIter;                           // itermax: maximum number of iteration
	int NonLinearMaxIter;                  // Maximum nonlinear iterations: number of multicorrection
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
	int *eqrho;
	NodeType *Node;
	ElementType *Element;
	double *F; // right hand side vector with size of the solution plus RANK before and after
	double *U; // vector of the solution with size of the solution plus RANK before and after
	double *u; // vector of the solution (current RANK)
	double *du; // Vector of differential (current RANK)
	double *f; // right hand side vector (current RANK)
	double *uB;
	double *duB;
	double *delta_old;
	AuxBuildStructuresType *AuxBuild;
	double *SendBuffer_bef[2], *SendBuffer_aft[2], *RecvBuffer_bef[2], *RecvBuffer_aft[2]; //Buffers to send and receive data to RANK-1 and RANK+1
	int *IdSend_bef, *IdSend_aft, *IdRecv_bef, *IdRecv_aft;//Vector to identify positions Node[I].invP_id in the buffers
}FemStructsType;

typedef struct
{
	double (*gamma)(double, double);
	double (*cv)(double, double);
	double (*rhopresc)(double, double);
	double (*v1presc)(double, double);
	double (*v2presc)(double, double);
	double (*epresc)(double, double);
	int (*InitialSolution)(ParametersType *, NodeType *, double *);

	double (*ShockCapture)(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4],
					double *, double, double, double, double, double, double, double, int, double *, double*);
	void (*Ax_Ay_calculations)(double, double, double [4], double [4][4], double [4][4]);
	double (*BC_theta)(double, double);
	void (*BC_no_penetrability)(int , int, int, NodeType *, double, double, double, double, double, double,
					    double, double, double, double, double [4][4], double [4][4], double [12][12], double [12][4],
						double [4][12], double [12], double [4], double [12], double [12], double [4], double [4]);
	void (*BC_General_theta)(int, int, int, NodeType *, double [3], double (*)(double, double));
	int (*StopCriteria)(ParametersType *, double, double, int);
	int (*StopTimeIntegration)(ParametersType *, double *, double *, double);

	void (*assembly)(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[12]);
	int (*mv)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precondR)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond_setup)(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
	int (*scaling)(ParametersType *, MatrixDataType *, FemStructsType *);
	int (*unscaling)(ParametersType *, MatrixDataType *, FemStructsType *, double *);
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

int setScaling(ParametersType *, FemFunctionsType *);

int setStabilizationForm(ParametersType *,FemFunctionsType *, FemOtherFunctionsType *, int (**)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *));

int setzeros(ParametersType *, MatrixDataType *);

double Delta_CAU(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4],
				double *, double, double, double, double, double, double, double, int, double *, double*);
				
double Delta_YZBeta(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4], 
				double *, double, double, double, double, double, double, double, int, double *, double *);

double Delta_YZBetaNMV(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4], 
				double *, double, double, double, double, double, double, double, int, double *, double *);

double Delta_DD(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4],
				double *, double, double, double, double, double, double, double, int, double *, double *);

void NO_BC_no_penetrability(int, int, int, NodeType *, double, double, double, double, double, double, double, double, double, double, 
				    double [4][4], double [4][4], double [12][12], double [12][4], double [4][12],
		         		double [12], double [4], double [12], double [12], double [4], double [4]);

void no_penetrability(double, double, double, double, double, double, double, double, double, double, double, double, 
	                 double, double, double, double, double [4][4], double [4][4], double [12][12], double [12][4], double [4][12],
		         double [12], double [4], double [12], double [12], double [4], double [4]);

void set_BC_no_penetrability(ParametersType *, FemFunctionsType *);

void csr_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[12]);

void csr_Initialization(ParametersType *, MatrixDataType *, FemStructsType *);

void csr_List_insertA(NodeListType **, int , int , int *);

int csr_search(int, int, NodeListType *);

void ebe_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[12]);

int Build_M_K_F_SUPG(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_M_K_F_DD_Transiente(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_M_F_DD_Transiente(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

void eval_U_dU(ParametersType *,FemStructsType *, FemFunctionsType *, double *,double *);

void dimensionless_Ax_Ay_calculations(double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4]);

void dimensional_Ax_Ay_calculations(double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4]);

void setDimensionlessness(ParametersType *, FemFunctionsType *);

int Predictor_Old(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_Old_BDF(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New_BDF(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New_TRBDF2(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_Old_TRBDF2(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Paraview_Output(int, int, ParametersType *, FemStructsType *, FemFunctionsType *);

int Paraview_Output_3D(int, int, ParametersType *, FemStructsType *, FemFunctionsType *);

int calculate_DaB(ParametersType *, FemStructsType *, double *, double *);

int uB_InitialSolution(ParametersType *, FemStructsType *, FemFunctionsType *, double *, double *);

void print_rho_Residue(FILE *outFile, int neqrho, double t,double *R, double *R_rho, int *eqrho);

void BC_theta_OK(int, int, int, NodeType *, double [3], double (*BC_theta)(double, double));

void BC_theta_NO(int, int, int, NodeType *, double [3], double (*BC_theta)(double, double));

int COMPARE_array_in_Preprocess (const void * a, const void * b);

void Permutation_of_LM(ParametersType *, MatrixDataType *, FemStructsType *);

void MPI_VectorUpdate(int, int, ParametersType *, FemStructsType *, double *, int);

#endif

