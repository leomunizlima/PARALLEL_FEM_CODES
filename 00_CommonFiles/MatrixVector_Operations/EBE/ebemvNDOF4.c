#include "../matvec.h"
#include "../../Allocation_Operations/allocations.h"

int ebemvNDOF4(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *P, double *Q)
{

	int I, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12;
	double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12;
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
		i10 = lm[I][9];
		i11 = lm[I][10];
		i12 = lm[I][11];

		p1 = P[i1];
		p2 = P[i2];
		p3 = P[i3];
		p4 = P[i4];
		p5 = P[i5];
		p6 = P[i6];
		p7 = P[i7];
		p8 = P[i8];
		p9 = P[i9];
		p10 = P[i10];
		p11 = P[i11];
		p12 = P[i12];

		Q[i1] += 	A[I][0]*p1 + A[I][1]*p2 + A[I][2]*p3
				+	A[I][3]*p4 + A[I][4]*p5 + A[I][5]*p6
				+	A[I][6]*p7 + A[I][7]*p8 + A[I][8]*p9
				+	A[I][9]*p10 + A[I][10]*p11 + A[I][11]*p12;

		Q[i2] += 	A[I][12]*p1 + A[I][13]*p2 + A[I][14]*p3
				+	A[I][15]*p4 + A[I][16]*p5 + A[I][17]*p6
				+	A[I][18]*p7 + A[I][19]*p8 + A[I][20]*p9
				+	A[I][21]*p10 + A[I][22]*p11 + A[I][23]*p12;

		Q[i3] += 	A[I][24]*p1 + A[I][25]*p2 + A[I][26]*p3
				+	A[I][27]*p4 + A[I][28]*p5 + A[I][29]*p6
				+	A[I][30]*p7 + A[I][31]*p8 + A[I][32]*p9
				+	A[I][33]*p10 + A[I][34]*p11 + A[I][35]*p12;

		Q[i4] += 	A[I][36]*p1 + A[I][37]*p2 + A[I][38]*p3
				+	A[I][39]*p4 + A[I][40]*p5 + A[I][41]*p6
				+	A[I][42]*p7 + A[I][43]*p8 + A[I][44]*p9
				+	A[I][45]*p10 + A[I][46]*p11 + A[I][47]*p12;

		Q[i5] += 	A[I][48]*p1 + A[I][49]*p2 + A[I][50]*p3
				+	A[I][51]*p4 + A[I][52]*p5 + A[I][53]*p6
				+	A[I][54]*p7 + A[I][55]*p8 + A[I][56]*p9
				+	A[I][57]*p10 + A[I][58]*p11 + A[I][59]*p12;

		Q[i6] += 	A[I][60]*p1 + A[I][61]*p2 + A[I][62]*p3
				+	A[I][63]*p4 + A[I][64]*p5 + A[I][65]*p6
				+	A[I][66]*p7 + A[I][67]*p8 + A[I][68]*p9
				+	A[I][69]*p10 + A[I][70]*p11 + A[I][71]*p12;

		Q[i7] += 	A[I][72]*p1 + A[I][73]*p2 + A[I][74]*p3
				+	A[I][75]*p4 + A[I][76]*p5 + A[I][77]*p6
				+	A[I][78]*p7 + A[I][79]*p8 + A[I][80]*p9
				+	A[I][81]*p10 + A[I][82]*p11 + A[I][83]*p12;

		Q[i8] += 	A[I][84]*p1 + A[I][85]*p2 + A[I][86]*p3
				+	A[I][87]*p4 + A[I][88]*p5 + A[I][89]*p6
				+	A[I][90]*p7 + A[I][91]*p8 + A[I][92]*p9
				+	A[I][93]*p10 + A[I][94]*p11 + A[I][95]*p12;

		Q[i9] += 	A[I][96]*p1 + A[I][97]*p2 + A[I][98]*p3
				+	A[I][99]*p4 + A[I][100]*p5 + A[I][101]*p6
				+	A[I][102]*p7 + A[I][103]*p8 + A[I][104]*p9
				+	A[I][105]*p10 + A[I][106]*p11 + A[I][107]*p12;

		Q[i10] += 	A[I][108]*p1 + A[I][109]*p2 + A[I][110]*p3
				+	A[I][111]*p4 + A[I][112]*p5 + A[I][113]*p6
				+	A[I][114]*p7 + A[I][115]*p8 + A[I][116]*p9
				+	A[I][117]*p10 + A[I][118]*p11 + A[I][119]*p12;

		Q[i11] += 	A[I][120]*p1 + A[I][121]*p2 + A[I][122]*p3
				+	A[I][123]*p4 + A[I][124]*p5 + A[I][125]*p6
				+	A[I][126]*p7 + A[I][127]*p8 + A[I][128]*p9
				+	A[I][129]*p10 + A[I][130]*p11 + A[I][131]*p12;

		Q[i12] += 	A[I][132]*p1 + A[I][133]*p2 + A[I][134]*p3
				+	A[I][135]*p4 + A[I][136]*p5 + A[I][137]*p6
				+	A[I][138]*p7 + A[I][139]*p8 + A[I][140]*p9
				+	A[I][141]*p10 + A[I][142]*p11 + A[I][143]*p12;

		Q[NEQ] = 0.0;
	}


	return 0;
}






