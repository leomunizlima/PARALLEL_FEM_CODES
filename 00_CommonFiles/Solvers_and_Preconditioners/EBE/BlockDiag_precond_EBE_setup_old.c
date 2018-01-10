#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int matrix_Inverse(double *, double *);
int fill_M(double [4][4],int , int , double *, int [12][12]);
int matrix_Bondary_Conditions(double [4][4], double [4][4], double [4][4], int, int);
int matrix_Preconditioning(double *, int, double *, double [4][4], double [4][4], double [4][4], int [12][12]);

int BlockBlockDiag_precond_EBE_setup(ParametersType *Parameters, MatrixDataType *MatrixData, NodeType *Node, ElementType *Element, int tag, double *F, int **lm)
{
	int I, J, E, I1, I2, I3;
	double **A = MatrixData->A[tag];
	int neq = Parameters-> neq;
	int nel = Parameters-> nel;
	int nnodes = Parameters->nnodes;
	int **Id, *IdAux;
	double **BlockDiag, *BlockDiagAux, **invBlockDiag, *invBlockDiagAux;

	int I_diag[3][48] = {{0,1,2,3,12,13,14,15,24,25,26,27,36,37,38,39},{52,53,54,55,64,65,66,67,76,77,78,79,88,89,90,91},{104,105,106,107,116,117,118,119,128,129,130,131,140,141,142,143}};

	int I_general[12][12] = {{0,1,2,3,4,5,6,7,8,9,10,11},{12,13,14,15,16,17,18,19,20,21,22,23},{24,25,26,27,28,29,30,31,32,33,34,35},{36,37,38,39,40,41,42,43,44,45,46,47},
				 {48,49,50,51,52,53,54,55,56,57,58,59},{60,61,62,63,64,65,66,67,68,69,70,71},{72,73,74,75,76,77,78,79,80,81,82,83},{84,85,86,87,88,89,90,91,92,93,94,95},
				 {96,97,98,99,100,101,102,103,104,105,106,107},{108,109,110,111,112,113,114,115,116,117,118,119},{120,121,122,123,124,125,126,127,128,129,130,131},
                                                                                                                                 {132,133,134,135,136,137,138,139,140,141,142,143}};	
	IdAux = (int*) mycalloc("IdAux of 'BlockBlockDiag_precond_EBE'",4*nnodes,sizeof(int));
	Id = (int**) mycalloc("Id of 'BlockBlockDiag_precond_EBE'", nnodes,sizeof(int*));
	BlockDiagAux = (double*) mycalloc("BlockDiagAux of 'BlockBlockDiag_precond_EBE'",16*nnodes,sizeof(double));
	BlockDiag = (double**) mycalloc("BlockDiag of 'BlockBlockDiag_precond_EBE'",nnodes,sizeof(double*));
	invBlockDiagAux = (double*) mycalloc("invBlockDiagAux of 'BlockBlockDiag_precond_EBE'",16*nnodes,sizeof(double));
	invBlockDiag = (double**) mycalloc("invBlockDiag of 'BlockBlockDiag_precond_EBE'",nnodes,sizeof(double*));

	for (I=0;I<nnodes;I++){
		BlockDiag[I] = &BlockDiagAux[16*I]; 
		invBlockDiag[I] = &invBlockDiagAux[16*I]; 
		Id[I] = &IdAux[4*I];
	}

	for (I=0;I<nnodes;I++){
		if (Node[I].id[0] < 0) 
			Id[I][0] = neq;
		else
			Id[I][0] = Node[I].id[0];

		if (Node[I].id[1] < 0) 
			Id[I][1] = neq; 
		else
			Id[I][1] = Node[I].id[1];
		if (Node[I].id[2] < 0) 
			Id[I][2] = neq; 
		else
			Id[I][2] = Node[I].id[2];
		if (Node[I].id[3] < 0)
			Id[I][3] = neq; 
		else
			Id[I][3] = Node[I].id[3];
	}

	for (E=0;E<nel;E++){
		I1 = Element[E].Vertex[0];	
		I2 = Element[E].Vertex[1];	
		I3 = Element[E].Vertex[2];	
		for (I=0;I<4;I++){
			for(J=0;J<4;J++){
				BlockDiag[I1][4*I+J] += A[E][I_diag[0][4*I+J]];
				BlockDiag[I2][4*I+J] += A[E][I_diag[1][4*I+J]];
				BlockDiag[I3][4*I+J] += A[E][I_diag[2][4*I+J]];
			}
		}
	}

	for (I=0;I<nnodes;I++)
		matrix_Inverse(BlockDiag[I],invBlockDiag[I]);
	

	double M11[4][4], M12[4][4], M13[4][4], M21[4][4], M22[4][4], M23[4][4], M31[4][4], M32[4][4], M33[4][4];
	
	for (E=0;E<nel;E++){
		I1 = Element[E].Vertex[0];	
		I2 = Element[E].Vertex[1];	
		I3 = Element[E].Vertex[2];	

		fill_M(M11,0,0,A[E],I_general);
		fill_M(M12,0,1,A[E],I_general);
		fill_M(M13,0,2,A[E],I_general);
		fill_M(M21,1,0,A[E],I_general);
		fill_M(M22,1,1,A[E],I_general);
		fill_M(M23,1,2,A[E],I_general);
		fill_M(M31,2,0,A[E],I_general);
		fill_M(M32,2,1,A[E],I_general);
		fill_M(M33,2,2,A[E],I_general);

		/* Adjusting boundary conditions on matrices M of the elements*/
		matrix_Bondary_Conditions(M11, M12, M13, Node[I1].id[0], 0);
		matrix_Bondary_Conditions(M11, M12, M13, Node[I1].id[1], 1);
		matrix_Bondary_Conditions(M11, M12, M13, Node[I1].id[2], 2);
		matrix_Bondary_Conditions(M11, M12, M13, Node[I1].id[3], 3);
			
		matrix_Bondary_Conditions(M22, M21, M23, Node[I2].id[0], 0);
		matrix_Bondary_Conditions(M22, M21, M23, Node[I2].id[1], 1);
		matrix_Bondary_Conditions(M22, M21, M23, Node[I2].id[2], 2);
		matrix_Bondary_Conditions(M22, M21, M23, Node[I2].id[3], 3);

		matrix_Bondary_Conditions(M33, M31, M32, Node[I3].id[0], 0);
		matrix_Bondary_Conditions(M33, M31, M32, Node[I3].id[1], 1);
		matrix_Bondary_Conditions(M33, M31, M32, Node[I3].id[2], 2);
		matrix_Bondary_Conditions(M33, M31, M32, Node[I3].id[3], 3);

		/* Preconditioning of the Ae = Inv*Me */	
		matrix_Preconditioning(A[E], 0, invBlockDiag[I1], M11, M12, M13, I_general);
		matrix_Preconditioning(A[E], 1, invBlockDiag[I2], M21, M22, M23, I_general);
		matrix_Preconditioning(A[E], 2, invBlockDiag[I3], M31, M32, M33, I_general);

	}

	/* Preconditioning of F */
	double a,b,c,d;

	for (I=0;I<nnodes;I++){
		a = F[Id[I][0]];
		b = F[Id[I][1]];
		c = F[Id[I][2]];
		d = F[Id[I][3]];
	
		F[Id[I][0]] = invBlockDiag[I][0]*a + invBlockDiag[I][1]*b + invBlockDiag[I][2]*c + invBlockDiag[I][3]*d;
		F[Id[I][1]] = invBlockDiag[I][4]*a + invBlockDiag[I][5]*b + invBlockDiag[I][6]*c + invBlockDiag[I][7]*d;
		F[Id[I][2]] = invBlockDiag[I][8]*a + invBlockDiag[I][9]*b + invBlockDiag[I][10]*c + invBlockDiag[I][11]*d;
		F[Id[I][3]] = invBlockDiag[I][12]*a + invBlockDiag[I][13]*b + invBlockDiag[I][14]*c + invBlockDiag[I][15]*d;

		F[neq] = 0;
	}
				
	free(BlockDiagAux);
	free(invBlockDiagAux);
	free(IdAux);
	free(BlockDiag);
	free(invBlockDiag);
	free(Id);

	return 0;
}

int matrix_Inverse(double *BlockDiag, double *invBlockDiag)
{
	double a, b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,det;

	/*This function calculates the inverse of a 4 x 4 matrix 
	
	   |a b c d|	
	A =|e f g h| 	
	   |i j k l|	
	   |m n o p|	*/

	a = BlockDiag[0];
	b = BlockDiag[1];
	c = BlockDiag[2];
	d = BlockDiag[3];
	e = BlockDiag[4];
	f = BlockDiag[5];
	g = BlockDiag[6];
	h = BlockDiag[7];
	i = BlockDiag[8];
	j = BlockDiag[9];
	k = BlockDiag[10];
	l = BlockDiag[11];
	m = BlockDiag[12];
	n = BlockDiag[13];
	o = BlockDiag[14];
	p = BlockDiag[15];
	
	/* Determinant calculations */
	det = (a*f*k*p - a*f*l*o - a*g*j*p + a*g*l*n + a*h*j*o - a*h*k*n
	      -b*e*k*p  + b*e*l*o + b*g*i*p - b*g*l*m - b*h*i*o + b*h*k*m +
	      c*(e*j*p - e*l*n - f*i*p + f*l*m + h*i*n - h*j*m) +
	      d*(-e*j*o + e*k*n + f*i*o - f*k*m - g*i*n + g*j*m));

//	if (fabs(det)<1e-16) printf("OPS!!!! Determante NULO!!!\n");

	double invdet = 1.0/det;

	/* Inverse of matrix calculations */
	invBlockDiag[0] = invdet*(-h*k*n + g*l*n + h*j*o - f*l*o - g*j*p + f*k*p);		 	
	invBlockDiag[1] = invdet*(d*k*n - c*l*n - d*j*o + b*l*o + c*j*p - b*k*p);		 	
	invBlockDiag[2] = invdet*(-d*g*n + c*h*n + d*f*o - b*h*o - c*f*p + b*g*p);		 	
	invBlockDiag[3] = invdet*(d*g*j - c*h*j - d*f*k + b*h*k + c*f*l - b*g*l);		 	
	invBlockDiag[4] = invdet*(h*k*m - g*l*m - h*i*o + e*l*o + g*i*p - e*k*p);		 	
	invBlockDiag[5] = invdet*(-d*k*m + c*l*m + d*i*o - a*l*o - c*i*p + a*k*p);
	invBlockDiag[6] = invdet*(d*g*m - c*h*m - d*e*o + a*h*o + c*e*p - a*g*p);		 	
	invBlockDiag[7] = invdet*(-d*g*i + c*h*i + d*e*k - a*h*k - c*e*l + a*g*l);
	invBlockDiag[8] = invdet*(-h*j*m + f*l*m + h*i*n - e*l*n - f*i*p + e*j*p);
	invBlockDiag[9] = invdet*(d*j*m - b*l*m - d*i*n + a*l*n + b*i*p - a*j*p);		 	
	invBlockDiag[10] = invdet*(-d*f*m + b*h*m + d*e*n - a*h*n - b*e*p + a*f*p);
	invBlockDiag[11] = invdet*(d*f*i - b*h*i - d*e*j + a*h*j + b*e*l - a*f*l);		 	
	invBlockDiag[12] = invdet*(g*j*m - f*k*m - g*i*n + e*k*n + f*i*o - e*j*o);		 	
	invBlockDiag[13] = invdet*(-c*j*m + b*k*m + c*i*n - a*k*n - b*i*o + a*j*o);
	invBlockDiag[14] = invdet*(c*j*m - b*g*m - c*e*n + a*g*n + b*e*o - a*f*o);		 	
	invBlockDiag[15] = invdet*(-c*f*i + b*g*i + c*e*j - a*g*j - b*e*k + a*f*k);

	return 0;
}		

int fill_M(double M[4][4],int I, int J, double *A, int I_g[12][12])
{
	M[0][0] = A[I_g[4*I][4*J]];
	M[0][1] = A[I_g[4*I][4*J+1]];
	M[0][2] = A[I_g[4*I][4*J+2]];
	M[0][3] = A[I_g[4*I][4*J+3]];
	M[1][0] = A[I_g[4*I+1][4*J]];
	M[1][1] = A[I_g[4*I+1][4*J+1]];
	M[1][2] = A[I_g[4*I+1][4*J+2]];
	M[1][3] = A[I_g[4*I+1][4*J+3]];
	M[2][0] = A[I_g[4*I+2][4*J]];
	M[2][1] = A[I_g[4*I+2][4*J+1]];
	M[2][2] = A[I_g[4*I+2][4*J+2]];
	M[2][3] = A[I_g[4*I+2][4*J+3]];
	M[3][0] = A[I_g[4*I+3][4*J]];
	M[3][1] = A[I_g[4*I+3][4*J+1]];
	M[3][2] = A[I_g[4*I+3][4*J+2]];
	M[3][3] = A[I_g[4*I+3][4*J+3]];

	return 0;
}


int matrix_Bondary_Conditions(double M1[4][4], double M2[4][4], double M3[4][4], int NodeId, int tag)
{
	if (NodeId < 0){
		M1[tag][0] = 0;
		M1[tag][1] = 0;
		M1[tag][2] = 0;
		M1[tag][3] = 0;
		M1[0][tag] = 0;
		M1[1][tag] = 0;
		M1[2][tag] = 0;
		M1[3][tag] = 0;
		M1[tag][tag] = 1;
		M2[tag][0] = 0;
		M2[tag][1] = 0;
		M2[tag][2] = 0;
		M2[tag][3] = 0;
		M3[tag][0] = 0;
		M3[tag][1] = 0;
		M3[tag][2] = 0;
		M3[tag][3] = 0;
	}

	return 0;	
}


int matrix_Preconditioning(double *A, int tag, double *invBlockDiag, double M1[4][4], double M2[4][4], double M3[4][4], int I_g[12][12])
{

	A[I_g[4*tag][0]]     = invBlockDiag[0]*M1[0][0]  + invBlockDiag[1]*M1[1][0]  + invBlockDiag[2]*M1[2][0]  + invBlockDiag[3]*M1[3][0];
	A[I_g[4*tag][1]]     = invBlockDiag[0]*M1[0][1]  + invBlockDiag[1]*M1[1][1]  + invBlockDiag[2]*M1[2][1]  + invBlockDiag[3]*M1[3][1];
	A[I_g[4*tag][2]]     = invBlockDiag[0]*M1[0][2]  + invBlockDiag[1]*M1[1][2]  + invBlockDiag[2]*M1[2][2]  + invBlockDiag[3]*M1[3][2];
	A[I_g[4*tag][3]]     = invBlockDiag[0]*M1[0][3]  + invBlockDiag[1]*M1[1][3]  + invBlockDiag[2]*M1[2][3]  + invBlockDiag[3]*M1[3][3];

	A[I_g[4*tag][4]]     = invBlockDiag[0]*M2[0][0]  + invBlockDiag[1]*M2[1][0]  + invBlockDiag[2]*M2[2][0]  + invBlockDiag[3]*M2[3][0];
	A[I_g[4*tag][5]]     = invBlockDiag[0]*M2[0][1]  + invBlockDiag[1]*M2[1][1]  + invBlockDiag[2]*M2[2][1]  + invBlockDiag[3]*M2[3][1];
	A[I_g[4*tag][6]]     = invBlockDiag[0]*M2[0][2]  + invBlockDiag[1]*M2[1][2]  + invBlockDiag[2]*M2[2][2]  + invBlockDiag[3]*M2[3][2];
	A[I_g[4*tag][7]]     = invBlockDiag[0]*M2[0][3]  + invBlockDiag[1]*M2[1][3]  + invBlockDiag[2]*M2[2][3]  + invBlockDiag[3]*M2[3][3];

	A[I_g[4*tag][8]]     = invBlockDiag[0]*M3[0][0]  + invBlockDiag[1]*M3[1][0]  + invBlockDiag[2]*M3[2][0]  + invBlockDiag[3]*M3[3][0];
	A[I_g[4*tag][9]]     = invBlockDiag[0]*M3[0][1]  + invBlockDiag[1]*M3[1][1]  + invBlockDiag[2]*M3[2][1]  + invBlockDiag[3]*M3[3][1];
	A[I_g[4*tag][10]]    = invBlockDiag[0]*M3[0][2]  + invBlockDiag[1]*M3[1][2]  + invBlockDiag[2]*M3[2][2]  + invBlockDiag[3]*M3[3][2];
	A[I_g[4*tag][11]]    = invBlockDiag[0]*M3[0][3]  + invBlockDiag[1]*M3[1][3]  + invBlockDiag[2]*M3[2][3]  + invBlockDiag[3]*M3[3][3];

	A[I_g[4*tag+1][0]]     = invBlockDiag[4]*M1[0][0]  + invBlockDiag[5]*M1[1][0]  + invBlockDiag[6]*M1[2][0]  + invBlockDiag[7]*M1[3][0];
	A[I_g[4*tag+1][1]]     = invBlockDiag[4]*M1[0][1]  + invBlockDiag[5]*M1[1][1]  + invBlockDiag[6]*M1[2][1]  + invBlockDiag[7]*M1[3][1];
	A[I_g[4*tag+1][2]]     = invBlockDiag[4]*M1[0][2]  + invBlockDiag[5]*M1[1][2]  + invBlockDiag[6]*M1[2][2]  + invBlockDiag[7]*M1[3][2];
	A[I_g[4*tag+1][3]]     = invBlockDiag[4]*M1[0][3]  + invBlockDiag[5]*M1[1][3]  + invBlockDiag[6]*M1[2][3]  + invBlockDiag[7]*M1[3][3];

	A[I_g[4*tag+1][4]]     = invBlockDiag[4]*M2[0][0]  + invBlockDiag[5]*M2[1][0]  + invBlockDiag[6]*M2[2][0]  + invBlockDiag[7]*M2[3][0];
	A[I_g[4*tag+1][5]]     = invBlockDiag[4]*M2[0][1]  + invBlockDiag[5]*M2[1][1]  + invBlockDiag[6]*M2[2][1]  + invBlockDiag[7]*M2[3][1];
	A[I_g[4*tag+1][6]]     = invBlockDiag[4]*M2[0][2]  + invBlockDiag[5]*M2[1][2]  + invBlockDiag[6]*M2[2][2]  + invBlockDiag[7]*M2[3][2];
	A[I_g[4*tag+1][7]]     = invBlockDiag[4]*M2[0][3]  + invBlockDiag[5]*M2[1][3]  + invBlockDiag[6]*M2[2][3]  + invBlockDiag[7]*M2[3][3];

	A[I_g[4*tag+1][8]]     = invBlockDiag[4]*M3[0][0]  + invBlockDiag[5]*M3[1][0]  + invBlockDiag[6]*M3[2][0]  + invBlockDiag[7]*M3[3][0];
	A[I_g[4*tag+1][9]]     = invBlockDiag[4]*M3[0][1]  + invBlockDiag[5]*M3[1][1]  + invBlockDiag[6]*M3[2][1]  + invBlockDiag[7]*M3[3][1];
	A[I_g[4*tag+1][10]]    = invBlockDiag[4]*M3[0][2]  + invBlockDiag[5]*M3[1][2]  + invBlockDiag[6]*M3[2][2]  + invBlockDiag[7]*M3[3][2];
	A[I_g[4*tag+1][11]]    = invBlockDiag[4]*M3[0][3]  + invBlockDiag[5]*M3[1][3]  + invBlockDiag[6]*M3[2][3]  + invBlockDiag[7]*M3[3][3];

	A[I_g[4*tag+2][0]]     = invBlockDiag[8]*M1[0][0]  + invBlockDiag[9]*M1[1][0]  + invBlockDiag[10]*M1[2][0]  + invBlockDiag[11]*M1[3][0];
	A[I_g[4*tag+2][1]]     = invBlockDiag[8]*M1[0][1]  + invBlockDiag[9]*M1[1][1]  + invBlockDiag[10]*M1[2][1]  + invBlockDiag[11]*M1[3][1];
	A[I_g[4*tag+2][2]]     = invBlockDiag[8]*M1[0][2]  + invBlockDiag[9]*M1[1][2]  + invBlockDiag[10]*M1[2][2]  + invBlockDiag[11]*M1[3][2];
	A[I_g[4*tag+2][3]]     = invBlockDiag[8]*M1[0][3]  + invBlockDiag[9]*M1[1][3]  + invBlockDiag[10]*M1[2][3]  + invBlockDiag[11]*M1[3][3];

	A[I_g[4*tag+2][4]]     = invBlockDiag[8]*M2[0][0]  + invBlockDiag[9]*M2[1][0]  + invBlockDiag[10]*M2[2][0]  + invBlockDiag[11]*M2[3][0];
	A[I_g[4*tag+2][5]]     = invBlockDiag[8]*M2[0][1]  + invBlockDiag[9]*M2[1][1]  + invBlockDiag[10]*M2[2][1]  + invBlockDiag[11]*M2[3][1];
	A[I_g[4*tag+2][6]]     = invBlockDiag[8]*M2[0][2]  + invBlockDiag[9]*M2[1][2]  + invBlockDiag[10]*M2[2][2]  + invBlockDiag[11]*M2[3][2];
	A[I_g[4*tag+2][7]]     = invBlockDiag[8]*M2[0][3]  + invBlockDiag[9]*M2[1][3]  + invBlockDiag[10]*M2[2][3]  + invBlockDiag[11]*M2[3][3];

	A[I_g[4*tag+2][8]]     = invBlockDiag[8]*M3[0][0]  + invBlockDiag[9]*M3[1][0]  + invBlockDiag[10]*M3[2][0]  + invBlockDiag[11]*M3[3][0];
	A[I_g[4*tag+2][9]]     = invBlockDiag[8]*M3[0][1]  + invBlockDiag[9]*M3[1][1]  + invBlockDiag[10]*M3[2][1]  + invBlockDiag[11]*M3[3][1];
	A[I_g[4*tag+2][10]]    = invBlockDiag[8]*M3[0][2]  + invBlockDiag[9]*M3[1][2]  + invBlockDiag[10]*M3[2][2]  + invBlockDiag[11]*M3[3][2];
	A[I_g[4*tag+2][11]]    = invBlockDiag[8]*M3[0][3]  + invBlockDiag[9]*M3[1][3]  + invBlockDiag[10]*M3[2][3]  + invBlockDiag[11]*M3[3][3];

	A[I_g[4*tag+3][0]]     = invBlockDiag[12]*M1[0][0]  + invBlockDiag[13]*M1[1][0]  + invBlockDiag[14]*M1[2][0]  + invBlockDiag[15]*M1[3][0];
	A[I_g[4*tag+3][1]]     = invBlockDiag[12]*M1[0][1]  + invBlockDiag[13]*M1[1][1]  + invBlockDiag[14]*M1[2][1]  + invBlockDiag[15]*M1[3][1];
	A[I_g[4*tag+3][2]]     = invBlockDiag[12]*M1[0][2]  + invBlockDiag[13]*M1[1][2]  + invBlockDiag[14]*M1[2][2]  + invBlockDiag[15]*M1[3][2];
	A[I_g[4*tag+3][3]]     = invBlockDiag[12]*M1[0][3]  + invBlockDiag[13]*M1[1][3]  + invBlockDiag[14]*M1[2][3]  + invBlockDiag[15]*M1[3][3];

	A[I_g[4*tag+3][4]]     = invBlockDiag[12]*M2[0][0]  + invBlockDiag[13]*M2[1][0]  + invBlockDiag[14]*M2[2][0]  + invBlockDiag[15]*M2[3][0];
	A[I_g[4*tag+3][5]]     = invBlockDiag[12]*M2[0][1]  + invBlockDiag[13]*M2[1][1]  + invBlockDiag[14]*M2[2][1]  + invBlockDiag[15]*M2[3][1];
	A[I_g[4*tag+3][6]]     = invBlockDiag[12]*M2[0][2]  + invBlockDiag[13]*M2[1][2]  + invBlockDiag[14]*M2[2][2]  + invBlockDiag[15]*M2[3][2];
	A[I_g[4*tag+3][7]]     = invBlockDiag[12]*M2[0][3]  + invBlockDiag[13]*M2[1][3]  + invBlockDiag[14]*M2[2][3]  + invBlockDiag[15]*M2[3][3];

	A[I_g[4*tag+3][8]]     = invBlockDiag[12]*M3[0][0]  + invBlockDiag[13]*M3[1][0]  + invBlockDiag[14]*M3[2][0]  + invBlockDiag[15]*M3[3][0];
	A[I_g[4*tag+3][9]]     = invBlockDiag[12]*M3[0][1]  + invBlockDiag[13]*M3[1][1]  + invBlockDiag[14]*M3[2][1]  + invBlockDiag[15]*M3[3][1];
	A[I_g[4*tag+3][10]]    = invBlockDiag[12]*M3[0][2]  + invBlockDiag[13]*M3[1][2]  + invBlockDiag[14]*M3[2][2]  + invBlockDiag[15]*M3[3][2];
	A[I_g[4*tag+3][11]]    = invBlockDiag[12]*M3[0][3]  + invBlockDiag[13]*M3[1][3]  + invBlockDiag[14]*M3[2][3]  + invBlockDiag[15]*M3[3][3];

	return 0;

}	






