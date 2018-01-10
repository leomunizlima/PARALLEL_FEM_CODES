#include "../matvec.h"

int ebemvNDOF3(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *P, double *Q)
{

	int I, i1, i2, i3, i4, i5, i6, i7, i8, i9;
	double p1, p2, p3, p4, p5, p6, p7, p8, p9;
	int NEQ, neq_bef, nel, RANK, NPROC;
	double **A;
	int **lm;

	NEQ = Parameters->NEQ;
	nel = Parameters->nel;
	neq_bef = Parameters->neq_bef;
	lm = FemStructs->lm;

	A = MatrixData->A;
	RANK = Parameters->RANK;
	NPROC = Parameters->NPROC;

	MPI_VectorUpdate(RANK,NPROC, Parameters, FemStructs, &P[neq_bef],0);

	for (I=0; I < NEQ; I++)
		Q[I] = 0;

	Q[NEQ] = 0;
	P[NEQ] = 0;

	for(I = 0; I < nel; I++){
		i1 = lm[I][0];
		i2 = lm[I][1];
		i3 = lm[I][2];
		i4 = lm[I][3];
		i5 = lm[I][4];
		i6 = lm[I][5];
		i7 = lm[I][6];
		i8 = lm[I][7];
		i9 = lm[I][8];
		
		p1 = P[i1];
		p2 = P[i2];
		p3 = P[i3];
		p4 = P[i4];
		p5 = P[i5];
		p6 = P[i6];
		p7 = P[i7];
		p8 = P[i8];
		p9 = P[i9];
		
		Q[i1] += 	A[I][0]*p1 + A[I][1]*p2 + A[I][2]*p3
				+	A[I][3]*p4 + A[I][4]*p5 + A[I][5]*p6
				+	A[I][6]*p7 + A[I][7]*p8 + A[I][8]*p9;
				

		Q[i2] += 	A[I][9]*p1 + A[I][10]*p2 + A[I][11]*p3
				+	A[I][12]*p4 + A[I][13]*p5 + A[I][14]*p6
				+	A[I][15]*p7 + A[I][16]*p8 + A[I][17]*p9;

		Q[i3] += 	A[I][18]*p1 + A[I][19]*p2 + A[I][20]*p3
				+	A[I][21]*p4 + A[I][22]*p5 + A[I][23]*p6
				+	A[I][24]*p7 + A[I][25]*p8 + A[I][26]*p9;

		Q[i4] += 	A[I][27]*p1 + A[I][28]*p2 + A[I][29]*p3
				+	A[I][30]*p4 + A[I][31]*p5 + A[I][32]*p6
				+	A[I][33]*p7 + A[I][34]*p8 + A[I][35]*p9;

		Q[i5] += 	A[I][36]*p1 + A[I][37]*p2 + A[I][38]*p3
				+	A[I][39]*p4 + A[I][40]*p5 + A[I][41]*p6
				+	A[I][42]*p7 + A[I][43]*p8 + A[I][44]*p9;

		Q[i6] += 	A[I][45]*p1 + A[I][46]*p2 + A[I][47]*p3
				+	A[I][48]*p4 + A[I][49]*p5 + A[I][50]*p6
				+	A[I][51]*p7 + A[I][52]*p8 + A[I][53]*p9;

		Q[i7] += 	A[I][54]*p1 + A[I][55]*p2 + A[I][56]*p3
				+	A[I][57]*p4 + A[I][58]*p5 + A[I][59]*p6
				+	A[I][60]*p7 + A[I][61]*p8 + A[I][62]*p9;

		Q[i8] += 	A[I][63]*p1 + A[I][64]*p2 + A[I][65]*p3
				+	A[I][66]*p4 + A[I][67]*p5 + A[I][68]*p6
				+	A[I][69]*p7 + A[I][70]*p8 + A[I][71]*p9;

		Q[i9] += 	A[I][72]*p1 + A[I][73]*p2 + A[I][74]*p3
				+	A[I][75]*p4 + A[I][76]*p5 + A[I][77]*p6
				+	A[I][78]*p7 + A[I][79]*p8 + A[I][80]*p9;

		Q[NEQ] = 0.0;

	}

	return 0;
}






