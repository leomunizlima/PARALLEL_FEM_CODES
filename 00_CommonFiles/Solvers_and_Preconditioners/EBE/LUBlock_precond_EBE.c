#include "../preconditioners.h"
/* Define preconditioner action on matrix-vector product */
/* p = A z
 * A = (LU)^{-1}
 * L * p = w
 * U * w = z
 * z out
*/
int LUBlock_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z){
	int I;
	int lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11; //auxiliar
	double z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11; //auxiliar
	int nel = Parameters->nel;  //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **LUe = MatrixData->LUe;

	z[neq] = 0;

	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];
		lm3 = lm[I][3];
		lm4 = lm[I][4];
		lm5 = lm[I][5];
		lm6 = lm[I][6];
		lm7 = lm[I][7];
		lm8 = lm[I][8];
		lm9 = lm[I][9];
		lm10 = lm[I][10];
		lm11 = lm[I][11];

		z0 = z[lm0];
		z1 = z[lm1];
		z2 = z[lm2];
		z3 = z[lm3];
		z4 = z[lm4];
		z5 = z[lm5];
		z6 = z[lm6];
		z7 = z[lm7];
		z8 = z[lm8];
		z9 = z[lm9];
		z10 = z[lm10];
		z11 = z[lm11];

		// L * p = w

		z4 += - LUe[I][48]*z0 - LUe[I][49]*z1 - LUe[I][50]*z2 - LUe[I][51]*z3; 
		z5 += - LUe[I][60]*z0 - LUe[I][61]*z1 - LUe[I][62]*z2 - LUe[I][63]*z3 
			- LUe[I][64]*z4;  
		z6 += - LUe[I][72]*z0 - LUe[I][73]*z1 - LUe[I][74]*z2 - LUe[I][75]*z3 
			- LUe[I][76]*z4 - LUe[I][77]*z5; 
		z7 += - LUe[I][84]*z0 - LUe[I][85]*z1 - LUe[I][86]*z2 - LUe[I][87]*z3 
			- LUe[I][88]*z4 - LUe[I][89]*z5 - LUe[I][90]*z6; 
		z8  += - LUe[I][96 ]*z0 - LUe[I][97 ]*z1 - LUe[I][98 ]*z2 - LUe[I][99 ]*z3 
		       - LUe[I][100]*z4 - LUe[I][101]*z5 - LUe[I][102]*z6 - LUe[I][103]*z7;
		z9  += - LUe[I][108]*z0 - LUe[I][109]*z1 - LUe[I][110]*z2 - LUe[I][111]*z3 
		       - LUe[I][112]*z4 - LUe[I][113]*z5 - LUe[I][114]*z6 - LUe[I][115]*z7 
		       - LUe[I][116]*z8;
		z10 += - LUe[I][120]*z0 - LUe[I][121]*z1 - LUe[I][122]*z2 - LUe[I][123]*z3 
		       - LUe[I][124]*z4 - LUe[I][125]*z5 - LUe[I][126]*z6 - LUe[I][127]*z7 
                       - LUe[I][128]*z8 - LUe[I][129]*z9;
		z11 += - LUe[I][132]*z0 - LUe[I][133]*z1 - LUe[I][134]*z2 - LUe[I][135]*z3 
		       - LUe[I][136]*z4 - LUe[I][137]*z5 - LUe[I][138]*z6 - LUe[I][139]*z7
                       - LUe[I][140]*z8 - LUe[I][141]*z9 - LUe[I][142]*z10;


		// U * w = z
		z11 = z11 * LUe[I][143];

		z10 = (z10 - LUe[I][131]*z11)*LUe[I][130];

		z9  = (z9  - LUe[I][118]*z10 - LUe[I][119]*z10)*LUe[I][117];

		z8  = (z8  - LUe[I][105]*z9 - LUe[I][106]*z10 - LUe[I][107]*z11)*LUe[I][104];

		z7  = (z7 - LUe[I][92]*z8 - LUe[I][93]*z9 - LUe[I][94]*z10 - LUe[I][95]*z11)*LUe[I][91];
		
		z6  = (z6 - LUe[I][79]*z7 - LUe[I][80]*z8 - LUe[I][81]*z9 - LUe[I][82]*z10 - LUe[I][83]*z11)*LUe[I][78];		
		z5  = (z5 - LUe[I][66]*z6 - LUe[I][67]*z7 - LUe[I][68]*z8 - LUe[I][69]*z9 - LUe[I][70]*z10 - 
		       LUe[I][71]*z11)*LUe[I][65];	

		z4  = (z4 - LUe[I][53]*z5 - LUe[I][54]*z6 - LUe[I][55]*z7 - LUe[I][56]*z8 - LUe[I][57]*z9 - 
		       LUe[I][58]*z10 - LUe[I][59]*z11)*LUe[I][52];		

		z3  = (z3 - LUe[I][40]*z4 - LUe[I][41]*z5 - LUe[I][42]*z6 - LUe[I][43]*z7 - LUe[I][44]*z8 - 
		       LUe[I][45]*z9 - LUe[I][46]*z10 - LUe[I][47]*z11)*LUe[I][39];	

		z2  = (z2 - LUe[I][27]*z3 - LUe[I][28]*z4 - LUe[I][29]*z5 - LUe[I][30]*z6 - LUe[I][31]*z7 - 
		       LUe[I][32]*z8 - LUe[I][33]*z9 - LUe[I][34]*z10 - LUe[I][35]*z11)*LUe[I][26];	

		z1 = (z1 - LUe[I][14]*z2 - LUe[I][15]*z3 - LUe[I][16]*z4 - LUe[I][17]*z5 - LUe[I][18]*z6 - 
		       LUe[I][19]*z7 - LUe[I][20]*z8 - LUe[I][21]*z9 - LUe[I][22]*z10 - LUe[I][23]*z11)*LUe[I][13];	

		z0 = (z0 - LUe[I][1]*z1 - LUe[I][2]*z2 - LUe[I][3]*z3 - LUe[I][4]*z4 - LUe[I][5]*z5 - LUe[I][6]*z6 - 
			LUe[I][7]*z7 - LUe[I][8]*z8 - LUe[I][9]*z9 - LUe[I][10]*z10 - LUe[I][11]*z11)*LUe[I][0];	
		

		z[lm0] = z0;
		z[lm1] = z1;
		z[lm2] = z2;
		z[lm3] = z3;
		z[lm4] = z4;
		z[lm5] = z5;
		z[lm6] = z6;
		z[lm7] = z7;
		z[lm8] = z8;
		z[lm9] = z9;
		z[lm10] = z10;
		z[lm11] = z11;

		z[neq] = 0.0;
	}

	return 0;
}


