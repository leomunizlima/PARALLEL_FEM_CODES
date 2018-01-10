#include "EulerEquations.h"

void no_penetrability(double s1, double s2, double s3, double c1, double c2, double c3, double alpha_dt, double delta, double deltaNMV, double Area, double y23, double y12, 
	                 double y31, double x32, double x21, double x13, double Ax[4][4], double Ay[4][4], double M1[12][12], double N1[12][4], double M2[4][12],
		         double Res1[12], double Res2[4], double Ue[12], double dUe[12], double UeB[4], double dUeB[4])

{
	//common parts
	double one_sixth = 1./6.;
	double nine_over_forty = 9./40.;
	double Area_over_six = Area/6.;
	double Area_over_twelve = Area/12.;
	double Area_3_over_twenty = 3.*Area/20.;
	double delta_quarter_Area = 0.25*delta*Area;
	double R1AxR1[4][4], R1AxR2[4][4], R1AxR3[4][4];	
	double R2AxR1[4][4], R2AxR2[4][4], R2AxR3[4][4];	
	double R3AxR1[4][4], R3AxR2[4][4], R3AxR3[4][4];	
	double R1AyR1[4][4], R1AyR2[4][4], R1AyR3[4][4];	
	double R2AyR1[4][4], R2AyR2[4][4], R2AyR3[4][4];	
	double R3AyR1[4][4], R3AyR2[4][4], R3AyR3[4][4];	
	double R1Ax[4][4], R2Ax[4][4], R3Ax[4][4];
	double R1Ay[4][4], R2Ay[4][4], R3Ay[4][4];
	double AxR1[4][4], AxR2[4][4], AxR3[4][4];
	double AyR1[4][4], AyR2[4][4], AyR3[4][4];
	double R1R1[4][4], R1R2[4][4], R1R3[4][4];
	double R2R2[4][4], R2R3[4][4];
	double R3R3[4][4];
	double R1[4][4], R2[4][4], R3[4][4];

	//R1AxR1
	R1AxR1[0][0] = Ax[0][0];
	R1AxR1[0][1] = Ax[0][2]*s1 + Ax[0][1]*c1;
	R1AxR1[0][2] = Ax[0][2]*c1 - Ax[0][1]*s1;
	R1AxR1[0][3] = Ax[0][3];
	R1AxR1[1][0] = Ax[1][0]*c1 - Ax[2][0]*s1;
	R1AxR1[1][1] = c1*(Ax[1][1]*c1 + Ax[1][2]*s1) - s1*(Ax[2][1]*c1 + Ax[2][2]*s1);
	R1AxR1[1][2] = c1*(Ax[1][2]*c1 - Ax[1][1]*s1) - s1*(Ax[2][2]*c1 - Ax[2][1]*s1);
	R1AxR1[1][3] = Ax[1][3]*c1 - Ax[2][3]*s1;
	R1AxR1[2][0] = Ax[1][0]*s1 + Ax[2][0]*c1;
	R1AxR1[2][1] = c1*(Ax[2][1]*c1 + Ax[2][2]*s1) + s1*(Ax[1][1]*c1 + Ax[1][2]*s1);
	R1AxR1[2][2] = c1*(Ax[2][2]*c1 - Ax[2][1]*s1) + s1*(Ax[1][2]*c1 - Ax[1][1]*s1);
	R1AxR1[2][3] = Ax[1][3]*s1 + Ax[2][3]*c1;
	R1AxR1[3][0] = Ax[3][0];
	R1AxR1[3][1] = Ax[3][2]*s1 + Ax[3][1]*c1;
	R1AxR1[3][2] = Ax[3][2]*c1 - Ax[3][1]*s1;
	R1AxR1[3][3] = Ax[3][3];

	//R1AxR2
	R1AxR2[0][0] = Ax[0][0];
	R1AxR2[0][1] = Ax[0][2]*s2 + Ax[0][1]*c2;
	R1AxR2[0][2] = Ax[0][2]*c2 - Ax[0][1]*s2;
	R1AxR2[0][3] = Ax[0][3];
	R1AxR2[1][0] = Ax[1][0]*c1 - Ax[2][0]*s1;
	R1AxR2[1][1] = c1*(Ax[1][1]*c2 + Ax[1][2]*s2) - s1*(Ax[2][1]*c2 + Ax[2][2]*s2);
	R1AxR2[1][2] = c1*(Ax[1][2]*c2 - Ax[1][1]*s2) - s1*(Ax[2][2]*c2 - Ax[2][1]*s2);
	R1AxR2[1][3] = Ax[1][3]*c1 - Ax[2][3]*s1;
	R1AxR2[2][0] = Ax[1][0]*s1 + Ax[2][0]*c1;
	R1AxR2[2][1] = c1*(Ax[2][1]*c2 + Ax[2][2]*s2) + s1*(Ax[1][1]*c2 + Ax[1][2]*s2);
	R1AxR2[2][2] = c1*(Ax[2][2]*c2 - Ax[2][1]*s2) + s1*(Ax[1][2]*c2 - Ax[1][1]*s2);
	R1AxR2[2][3] = Ax[1][3]*s1 + Ax[2][3]*c1;
	R1AxR2[3][0] = Ax[3][0];
	R1AxR2[3][1] = Ax[3][2]*s2 + Ax[3][1]*c2;
	R1AxR2[3][2] = Ax[3][2]*c2 - Ax[3][1]*s2;
	R1AxR2[3][3] = Ax[3][3];

	//R1AxR3
	R1AxR3[0][0] = Ax[0][0];
	R1AxR3[0][1] = Ax[0][2]*s3 + Ax[0][1]*c3;
	R1AxR3[0][2] = Ax[0][2]*c3 - Ax[0][1]*s3;
	R1AxR3[0][3] = Ax[0][3];
	R1AxR3[1][0] = Ax[1][0]*c1 - Ax[2][0]*s1;
	R1AxR3[1][1] = c1*(Ax[1][1]*c3 + Ax[1][2]*s3) - s1*(Ax[2][1]*c3 + Ax[2][2]*s3);
	R1AxR3[1][2] = c1*(Ax[1][2]*c3 - Ax[1][1]*s3) - s1*(Ax[2][2]*c3 - Ax[2][1]*s3);
	R1AxR3[1][3] = Ax[1][3]*c1 - Ax[2][3]*s1;
	R1AxR3[2][0] = Ax[1][0]*s1 + Ax[2][0]*c1;
	R1AxR3[2][1] = c1*(Ax[2][1]*c3 + Ax[2][2]*s3) + s1*(Ax[1][1]*c3 + Ax[1][2]*s3);
	R1AxR3[2][2] = c1*(Ax[2][2]*c3 - Ax[2][1]*s3) + s1*(Ax[1][2]*c3 - Ax[1][1]*s3);
	R1AxR3[2][3] = Ax[1][3]*s1 + Ax[2][3]*c1;
	R1AxR3[3][0] = Ax[3][0];
	R1AxR3[3][1] = Ax[3][2]*s3 + Ax[3][1]*c3;
	R1AxR3[3][2] = Ax[3][2]*c3 - Ax[3][1]*s3;
	R1AxR3[3][3] = Ax[3][3];

	//R2AxR1
	R2AxR1[0][0] = Ax[0][0];
	R2AxR1[0][1] = Ax[0][2]*s1 + Ax[0][1]*c1;
	R2AxR1[0][2] = Ax[0][2]*c1 - Ax[0][1]*s1;
	R2AxR1[0][3] = Ax[0][3];
	R2AxR1[1][0] = Ax[1][0]*c2 - Ax[2][0]*s2;
	R2AxR1[1][1] = c2*(Ax[1][1]*c1 + Ax[1][2]*s1) - s2*(Ax[2][1]*c1 + Ax[2][2]*s1);
	R2AxR1[1][2] = c2*(Ax[1][2]*c1 - Ax[1][1]*s1) - s2*(Ax[2][2]*c1 - Ax[2][1]*s1);
	R2AxR1[1][3] = Ax[1][3]*c2 - Ax[2][3]*s2;
	R2AxR1[2][0] = Ax[1][0]*s2 + Ax[2][0]*c2;
	R2AxR1[2][1] = c2*(Ax[2][1]*c1 + Ax[2][2]*s1) + s2*(Ax[1][1]*c1 + Ax[1][2]*s1);
	R2AxR1[2][2] = c2*(Ax[2][2]*c1 - Ax[2][1]*s1) + s2*(Ax[1][2]*c1 - Ax[1][1]*s1);
	R2AxR1[2][3] = Ax[1][3]*s2 + Ax[2][3]*c2;
	R2AxR1[3][0] = Ax[3][0];
	R2AxR1[3][1] = Ax[3][2]*s1 + Ax[3][1]*c1;
	R2AxR1[3][2] = Ax[3][2]*c1 - Ax[3][1]*s1;
	R2AxR1[3][3] = Ax[3][3];

	//R2AxR2
	R2AxR2[0][0] = Ax[0][0];
	R2AxR2[0][1] = Ax[0][2]*s2 + Ax[0][1]*c2;
	R2AxR2[0][2] = Ax[0][2]*c2 - Ax[0][1]*s2;
	R2AxR2[0][3] = Ax[0][3];
	R2AxR2[1][0] = Ax[1][0]*c2 - Ax[2][0]*s2;
	R2AxR2[1][1] = c2*(Ax[1][1]*c2 + Ax[1][2]*s2) - s2*(Ax[2][1]*c2 + Ax[2][2]*s2);
	R2AxR2[1][2] = c2*(Ax[1][2]*c2 - Ax[1][1]*s2) - s2*(Ax[2][2]*c2 - Ax[2][1]*s2);
	R2AxR2[1][3] = Ax[1][3]*c2 - Ax[2][3]*s2;
	R2AxR2[2][0] = Ax[1][0]*s2 + Ax[2][0]*c2;
	R2AxR2[2][1] = c2*(Ax[2][1]*c2 + Ax[2][2]*s2) + s2*(Ax[1][1]*c2 + Ax[1][2]*s2);
	R2AxR2[2][2] = c2*(Ax[2][2]*c2 - Ax[2][1]*s2) + s2*(Ax[1][2]*c2 - Ax[1][1]*s2);
	R2AxR2[2][3] = Ax[1][3]*s2 + Ax[2][3]*c2;
	R2AxR2[3][0] = Ax[3][0];
	R2AxR2[3][1] = Ax[3][2]*s2 + Ax[3][1]*c2;
	R2AxR2[3][2] = Ax[3][2]*c2 - Ax[3][1]*s2;
	R2AxR2[3][3] = Ax[3][3];

	//R2AxR3	
	R2AxR3[0][0] = Ax[0][0];
	R2AxR3[0][1] = Ax[0][2]*s3 + Ax[0][1]*c3;
	R2AxR3[0][2] = Ax[0][2]*c3 - Ax[0][1]*s3;
	R2AxR3[0][3] = Ax[0][3];
	R2AxR3[1][0] = Ax[1][0]*c2 - Ax[2][0]*s2;
	R2AxR3[1][1] = c2*(Ax[1][1]*c3 + Ax[1][2]*s3) - s2*(Ax[2][1]*c3 + Ax[2][2]*s3);
	R2AxR3[1][2] = c2*(Ax[1][2]*c3 - Ax[1][1]*s3) - s2*(Ax[2][2]*c3 - Ax[2][1]*s3);
	R2AxR3[1][3] = Ax[1][3]*c2 - Ax[2][3]*s2;
	R2AxR3[2][0] = Ax[1][0]*s2 + Ax[2][0]*c2;
	R2AxR3[2][1] = c2*(Ax[2][1]*c3 + Ax[2][2]*s3) + s2*(Ax[1][1]*c3 + Ax[1][2]*s3);
	R2AxR3[2][2] = c2*(Ax[2][2]*c3 - Ax[2][1]*s3) + s2*(Ax[1][2]*c3 - Ax[1][1]*s3);
	R2AxR3[2][3] = Ax[1][3]*s2 + Ax[2][3]*c2;
	R2AxR3[3][0] = Ax[3][0];
	R2AxR3[3][1] = Ax[3][2]*s3 + Ax[3][1]*c3;
	R2AxR3[3][2] = Ax[3][2]*c3 - Ax[3][1]*s3;
	R2AxR3[3][3] = Ax[3][3];

	//R3AxR1
	R3AxR1[0][0] = Ax[0][0];
	R3AxR1[0][1] = Ax[0][2]*s1 + Ax[0][1]*c1;
	R3AxR1[0][2] = Ax[0][2]*c1 - Ax[0][1]*s1;
	R3AxR1[0][3] = Ax[0][3];
	R3AxR1[1][0] = Ax[1][0]*c3 - Ax[2][0]*s3;
	R3AxR1[1][1] = c3*(Ax[1][1]*c1 + Ax[1][2]*s1) - s3*(Ax[2][1]*c1 + Ax[2][2]*s1);
	R3AxR1[1][2] = c3*(Ax[1][2]*c1 - Ax[1][1]*s1) - s3*(Ax[2][2]*c1 - Ax[2][1]*s1);
	R3AxR1[1][3] = Ax[1][3]*c3 - Ax[2][3]*s3;
	R3AxR1[2][0] = Ax[1][0]*s3 + Ax[2][0]*c3;
	R3AxR1[2][1] = c3*(Ax[2][1]*c1 + Ax[2][2]*s1) + s3*(Ax[1][1]*c1 + Ax[1][2]*s1);
	R3AxR1[2][2] = c3*(Ax[2][2]*c1 - Ax[2][1]*s1) + s3*(Ax[1][2]*c1 - Ax[1][1]*s1);
	R3AxR1[2][3] = Ax[1][3]*s3 + Ax[2][3]*c3;
	R3AxR1[3][0] = Ax[3][0];
	R3AxR1[3][1] = Ax[3][2]*s1 + Ax[3][1]*c1;
	R3AxR1[3][2] = Ax[3][2]*c1 - Ax[3][1]*s1;
	R3AxR1[3][3] = Ax[3][3];

	//R3AxR2	
	R3AxR2[0][0] = Ax[0][0];
	R3AxR2[0][1] = Ax[0][2]*s2 + Ax[0][1]*c2;
	R3AxR2[0][2] = Ax[0][2]*c2 - Ax[0][1]*s2;
	R3AxR2[0][3] = Ax[0][3];
	R3AxR2[1][0] = Ax[1][0]*c3 - Ax[2][0]*s3;
	R3AxR2[1][1] = c3*(Ax[1][1]*c2 + Ax[1][2]*s2) - s3*(Ax[2][1]*c2 + Ax[2][2]*s2);
	R3AxR2[1][2] = c3*(Ax[1][2]*c2 - Ax[1][1]*s2) - s3*(Ax[2][2]*c2 - Ax[2][1]*s2);
	R3AxR2[1][3] = Ax[1][3]*c3 - Ax[2][3]*s3;
	R3AxR2[2][0] = Ax[1][0]*s3 + Ax[2][0]*c3;
	R3AxR2[2][1] = c3*(Ax[2][1]*c2 + Ax[2][2]*s2) + s3*(Ax[1][1]*c2 + Ax[1][2]*s2);
	R3AxR2[2][2] = c3*(Ax[2][2]*c2 - Ax[2][1]*s2) + s3*(Ax[1][2]*c2 - Ax[1][1]*s2);
	R3AxR2[2][3] = Ax[1][3]*s3 + Ax[2][3]*c3;
	R3AxR2[3][0] = Ax[3][0];
	R3AxR2[3][1] = Ax[3][2]*s2 + Ax[3][1]*c2;
	R3AxR2[3][2] = Ax[3][2]*c2 - Ax[3][1]*s2;
	R3AxR2[3][3] = Ax[3][3];

	//R3AxR3	
	R3AxR3[0][0] = Ax[0][0];
	R3AxR3[0][1] = Ax[0][2]*s3 + Ax[0][1]*c3;
	R3AxR3[0][2] = Ax[0][2]*c3 - Ax[0][1]*s3;
	R3AxR3[0][3] = Ax[0][3];
	R3AxR3[1][0]= Ax[1][0]*c3 - Ax[2][0]*s3;
	R3AxR3[1][1]= c3*(Ax[1][1]*c3 + Ax[1][2]*s3) - s3*(Ax[2][1]*c3 + Ax[2][2]*s3);
	R3AxR3[1][2]= c3*(Ax[1][2]*c3 - Ax[1][1]*s3) - s3*(Ax[2][2]*c3 - Ax[2][1]*s3);
	R3AxR3[1][3]= Ax[1][3]*c3 - Ax[2][3]*s3;
	R3AxR3[2][0]= Ax[1][0]*s3 + Ax[2][0]*c3;
	R3AxR3[2][1]= c3*(Ax[2][1]*c3 + Ax[2][2]*s3) + s3*(Ax[1][1]*c3 + Ax[1][2]*s3);
	R3AxR3[2][2]= c3*(Ax[2][2]*c3 - Ax[2][1]*s3) + s3*(Ax[1][2]*c3 - Ax[1][1]*s3);
	R3AxR3[2][3]= Ax[1][3]*s3 + Ax[2][3]*c3;
	R3AxR3[3][0] = Ax[3][0];
	R3AxR3[3][1]= Ax[3][2]*s3 + Ax[3][1]*c3;
	R3AxR3[3][2]= Ax[3][2]*c3 - Ax[3][1]*s3;
	R3AxR3[3][3] = Ax[3][3];

		
	//R1Ax
	R1Ax[0][0] = Ax[0][0];
	R1Ax[0][1] = Ax[0][1];
	R1Ax[0][2] = Ax[0][2];
	R1Ax[0][3] = Ax[0][3];
	R1Ax[1][0] = Ax[1][0]*c1 - Ax[2][0]*s1;
	R1Ax[1][1] = Ax[1][1]*c1 - Ax[2][1]*s1;
	R1Ax[1][2] = Ax[1][2]*c1 - Ax[2][2]*s1;
	R1Ax[1][3] = Ax[1][3]*c1 - Ax[2][3]*s1;
	R1Ax[2][0] = Ax[1][0]*s1 + Ax[2][0]*c1;
	R1Ax[2][1] = Ax[1][1]*s1 + Ax[2][1]*c1;
	R1Ax[2][2] = Ax[1][2]*s1 + Ax[2][2]*c1;
	R1Ax[2][3] = Ax[1][3]*s1 + Ax[2][3]*c1;
	R1Ax[3][0] = Ax[3][0];
	R1Ax[3][1] = Ax[3][1];
	R1Ax[3][2] = Ax[3][2];
	R1Ax[3][3] = Ax[3][3];
	
	//R2Ax
	R2Ax[0][0] = Ax[0][0];
	R2Ax[0][1] = Ax[0][1];
	R2Ax[0][2] = Ax[0][2];
	R2Ax[0][3] = Ax[0][3];
	R2Ax[1][0] = Ax[1][0]*c2 - Ax[2][0]*s2;
	R2Ax[1][1] = Ax[1][1]*c2 - Ax[2][1]*s2;
	R2Ax[1][2] = Ax[1][2]*c2 - Ax[2][2]*s2;
	R2Ax[1][3] = Ax[1][3]*c2 - Ax[2][3]*s2;
	R2Ax[2][0] = Ax[1][0]*s2 + Ax[2][0]*c2;
	R2Ax[2][1] = Ax[1][1]*s2 + Ax[2][1]*c2;
	R2Ax[2][2] = Ax[1][2]*s2 + Ax[2][2]*c2;
	R2Ax[2][3] = Ax[1][3]*s2 + Ax[2][3]*c2;
	R2Ax[3][0] = Ax[3][0];
	R2Ax[3][1] = Ax[3][1];
	R2Ax[3][2] = Ax[3][2];
	R2Ax[3][3] = Ax[3][3];
	
	//R3Ax
	R3Ax[0][0] = Ax[0][0];
	R3Ax[0][1] = Ax[0][1];
	R3Ax[0][2] = Ax[0][2];
	R3Ax[0][3] = Ax[0][3];
	R3Ax[1][0] = Ax[1][0]*c3 - Ax[2][0]*s3;
	R3Ax[1][1] = Ax[1][1]*c3 - Ax[2][1]*s3;
	R3Ax[1][2] = Ax[1][2]*c3 - Ax[2][2]*s3;
	R3Ax[1][3] = Ax[1][3]*c3 - Ax[2][3]*s3;
	R3Ax[2][0] = Ax[1][0]*s3 + Ax[2][0]*c3;
	R3Ax[2][1] = Ax[1][1]*s3 + Ax[2][1]*c3;
	R3Ax[2][2] = Ax[1][2]*s3 + Ax[2][2]*c3;
	R3Ax[2][3] = Ax[1][3]*s3 + Ax[2][3]*c3;
	R3Ax[3][0] = Ax[3][0];
	R3Ax[3][1] = Ax[3][1];
	R3Ax[3][2] = Ax[3][2];
	R3Ax[3][3] = Ax[3][3];

	//AxR1
	AxR1[0][0] = Ax[0][0];
	AxR1[0][1] = Ax[0][2]*s1 + Ax[0][1]*c1;
	AxR1[0][2] = Ax[0][2]*c1 - Ax[0][1]*s1;
	AxR1[0][3] = Ax[0][3];	
	AxR1[1][0] = Ax[1][0];	
	AxR1[1][1] = Ax[1][2]*s1 + Ax[1][1]*c1;
	AxR1[1][2] = Ax[1][2]*c1 - Ax[1][1]*s1;
	AxR1[1][3] = Ax[1][3];	
	AxR1[2][0] = Ax[2][0];	
	AxR1[2][1] = Ax[2][2]*s1 + Ax[2][1]*c1;
	AxR1[2][2] = Ax[2][2]*c1 - Ax[2][1]*s1;
	AxR1[2][3] = Ax[2][3];	
	AxR1[3][0] = Ax[3][0];	
	AxR1[3][1] = Ax[3][2]*s1 + Ax[3][1]*c1;
	AxR1[3][2] = Ax[3][2]*c1 - Ax[3][1]*s1;
	AxR1[3][3] = Ax[3][3];	
	
	//AxR2
	AxR2[0][0] = Ax[0][0];
	AxR2[0][1] = Ax[0][2]*s2 + Ax[0][1]*c2;
	AxR2[0][2] = Ax[0][2]*c2 - Ax[0][1]*s2;
	AxR2[0][3] = Ax[0][3];	
	AxR2[1][0] = Ax[1][0];	
	AxR2[1][1] = Ax[1][2]*s2 + Ax[1][1]*c2;
	AxR2[1][2] = Ax[1][2]*c2 - Ax[1][1]*s2;
	AxR2[1][3] = Ax[1][3];	
	AxR2[2][0] = Ax[2][0];	
	AxR2[2][1] = Ax[2][2]*s2 + Ax[2][1]*c2;
	AxR2[2][2] = Ax[2][2]*c2 - Ax[2][1]*s2;
	AxR2[2][3] = Ax[2][3];	
	AxR2[3][0] = Ax[3][0];	
	AxR2[3][1] = Ax[3][2]*s2 + Ax[3][1]*c2;
	AxR2[3][2] = Ax[3][2]*c2 - Ax[3][1]*s2;
	AxR2[3][3] = Ax[3][3];	

	//AxR3
	AxR3[0][0] = Ax[0][0];
	AxR3[0][1] = Ax[0][2]*s3 + Ax[0][1]*c3;
	AxR3[0][2] = Ax[0][2]*c3 - Ax[0][1]*s3;
	AxR3[0][3] = Ax[0][3];	
	AxR3[1][0] = Ax[1][0];	
	AxR3[1][1] = Ax[1][2]*s3 + Ax[1][1]*c3;
	AxR3[1][2] = Ax[1][2]*c3 - Ax[1][1]*s3;
	AxR3[1][3] = Ax[1][3];	
	AxR3[2][0] = Ax[2][0];	
	AxR3[2][1] = Ax[2][2]*s3 + Ax[2][1]*c3;
	AxR3[2][2] = Ax[2][2]*c3 - Ax[2][1]*s3;
	AxR3[2][3] = Ax[2][3];	
	AxR3[3][0] = Ax[3][0];	
	AxR3[3][1] = Ax[3][2]*s3 + Ax[3][1]*c3;
	AxR3[3][2] = Ax[3][2]*c3 - Ax[3][1]*s3;
	AxR3[3][3] = Ax[3][3];	

	//R1AyR1
	R1AyR1[0][0] = Ay[0][0];
	R1AyR1[0][1] = Ay[0][2]*s1 + Ay[0][1]*c1;
	R1AyR1[0][2] = Ay[0][2]*c1 - Ay[0][1]*s1;
	R1AyR1[0][3] = Ay[0][3];
	R1AyR1[1][0] = Ay[1][0]*c1 - Ay[2][0]*s1;
	R1AyR1[1][1] = c1*(Ay[1][1]*c1 + Ay[1][2]*s1) - s1*(Ay[2][1]*c1 + Ay[2][2]*s1);
	R1AyR1[1][2] = c1*(Ay[1][2]*c1 - Ay[1][1]*s1) - s1*(Ay[2][2]*c1 - Ay[2][1]*s1);
	R1AyR1[1][3] = Ay[1][3]*c1 - Ay[2][3]*s1;
	R1AyR1[2][0] = Ay[1][0]*s1 + Ay[2][0]*c1;
	R1AyR1[2][1] = c1*(Ay[2][1]*c1 + Ay[2][2]*s1) + s1*(Ay[1][1]*c1 + Ay[1][2]*s1);
	R1AyR1[2][2] = c1*(Ay[2][2]*c1 - Ay[2][1]*s1) + s1*(Ay[1][2]*c1 - Ay[1][1]*s1);
	R1AyR1[2][3] = Ay[1][3]*s1 + Ay[2][3]*c1;
	R1AyR1[3][0] = Ay[3][0];
	R1AyR1[3][1] = Ay[3][2]*s1 + Ay[3][1]*c1;
	R1AyR1[3][2] = Ay[3][2]*c1 - Ay[3][1]*s1;
	R1AyR1[3][3] = Ay[3][3];

	//R1AyR2
	R1AyR2[0][0] = Ay[0][0];
	R1AyR2[0][1] = Ay[0][2]*s2 + Ay[0][1]*c2;
	R1AyR2[0][2] = Ay[0][2]*c2 - Ay[0][1]*s2;
	R1AyR2[0][3] = Ay[0][3];
	R1AyR2[1][0] = Ay[1][0]*c1 - Ay[2][0]*s1;
	R1AyR2[1][1] = c1*(Ay[1][1]*c2 + Ay[1][2]*s2) - s1*(Ay[2][1]*c2 + Ay[2][2]*s2);
	R1AyR2[1][2] = c1*(Ay[1][2]*c2 - Ay[1][1]*s2) - s1*(Ay[2][2]*c2 - Ay[2][1]*s2);
	R1AyR2[1][3] = Ay[1][3]*c1 - Ay[2][3]*s1;
	R1AyR2[2][0] = Ay[1][0]*s1 + Ay[2][0]*c1;
	R1AyR2[2][1] = c1*(Ay[2][1]*c2 + Ay[2][2]*s2) + s1*(Ay[1][1]*c2 + Ay[1][2]*s2);
	R1AyR2[2][2] = c1*(Ay[2][2]*c2 - Ay[2][1]*s2) + s1*(Ay[1][2]*c2 - Ay[1][1]*s2);
	R1AyR2[2][3] = Ay[1][3]*s1 + Ay[2][3]*c1;
	R1AyR2[3][0] = Ay[3][0];
	R1AyR2[3][1] = Ay[3][2]*s2 + Ay[3][1]*c2;
	R1AyR2[3][2] = Ay[3][2]*c2 - Ay[3][1]*s2;
	R1AyR2[3][3] = Ay[3][3];

	//R1AyR3
	R1AyR3[0][0] = Ay[0][0];
	R1AyR3[0][1] = Ay[0][2]*s3 + Ay[0][1]*c3;
	R1AyR3[0][2] = Ay[0][2]*c3 - Ay[0][1]*s3;
	R1AyR3[0][3] = Ay[0][3];
	R1AyR3[1][0] = Ay[1][0]*c1 - Ay[2][0]*s1;
	R1AyR3[1][1] = c1*(Ay[1][1]*c3 + Ay[1][2]*s3) - s1*(Ay[2][1]*c3 + Ay[2][2]*s3);
	R1AyR3[1][2] = c1*(Ay[1][2]*c3 - Ay[1][1]*s3) - s1*(Ay[2][2]*c3 - Ay[2][1]*s3);
	R1AyR3[1][3] = Ay[1][3]*c1 - Ay[2][3]*s1;
	R1AyR3[2][0] = Ay[1][0]*s1 + Ay[2][0]*c1;
	R1AyR3[2][1] = c1*(Ay[2][1]*c3 + Ay[2][2]*s3) + s1*(Ay[1][1]*c3 + Ay[1][2]*s3);
	R1AyR3[2][2] = c1*(Ay[2][2]*c3 - Ay[2][1]*s3) + s1*(Ay[1][2]*c3 - Ay[1][1]*s3);
	R1AyR3[2][3] = Ay[1][3]*s1 + Ay[2][3]*c1;
	R1AyR3[3][0] = Ay[3][0];
	R1AyR3[3][1] = Ay[3][2]*s3 + Ay[3][1]*c3;
	R1AyR3[3][2] = Ay[3][2]*c3 - Ay[3][1]*s3;
	R1AyR3[3][3] = Ay[3][3];

	//R2AyR1
	R2AyR1[0][0] = Ay[0][0];
	R2AyR1[0][1] = Ay[0][2]*s1 + Ay[0][1]*c1;
	R2AyR1[0][2] = Ay[0][2]*c1 - Ay[0][1]*s1;
	R2AyR1[0][3] = Ay[0][3];
	R2AyR1[1][0] = Ay[1][0]*c2 - Ay[2][0]*s2;
	R2AyR1[1][1] = c2*(Ay[1][1]*c1 + Ay[1][2]*s1) - s2*(Ay[2][1]*c1 + Ay[2][2]*s1);
	R2AyR1[1][2] = c2*(Ay[1][2]*c1 - Ay[1][1]*s1) - s2*(Ay[2][2]*c1 - Ay[2][1]*s1);
	R2AyR1[1][3] = Ay[1][3]*c2 - Ay[2][3]*s2;
	R2AyR1[2][0] = Ay[1][0]*s2 + Ay[2][0]*c2;
	R2AyR1[2][1] = c2*(Ay[2][1]*c1 + Ay[2][2]*s1) + s2*(Ay[1][1]*c1 + Ay[1][2]*s1);
	R2AyR1[2][2] = c2*(Ay[2][2]*c1 - Ay[2][1]*s1) + s2*(Ay[1][2]*c1 - Ay[1][1]*s1);
	R2AyR1[2][3] = Ay[1][3]*s2 + Ay[2][3]*c2;
	R2AyR1[3][0] = Ay[3][0];
	R2AyR1[3][1] = Ay[3][2]*s1 + Ay[3][1]*c1;
	R2AyR1[3][2] = Ay[3][2]*c1 - Ay[3][1]*s1;
	R2AyR1[3][3] = Ay[3][3];

	//R2AyR2
	R2AyR2[0][0] = Ay[0][0];
	R2AyR2[0][1] = Ay[0][2]*s2 + Ay[0][1]*c2;
	R2AyR2[0][2] = Ay[0][2]*c2 - Ay[0][1]*s2;
	R2AyR2[0][3] = Ay[0][3];
	R2AyR2[1][0] = Ay[1][0]*c2 - Ay[2][0]*s2;
	R2AyR2[1][1] = c2*(Ay[1][1]*c2 + Ay[1][2]*s2) - s2*(Ay[2][1]*c2 + Ay[2][2]*s2);
	R2AyR2[1][2] = c2*(Ay[1][2]*c2 - Ay[1][1]*s2) - s2*(Ay[2][2]*c2 - Ay[2][1]*s2);
	R2AyR2[1][3] = Ay[1][3]*c2 - Ay[2][3]*s2;
	R2AyR2[2][0] = Ay[1][0]*s2 + Ay[2][0]*c2;
	R2AyR2[2][1] = c2*(Ay[2][1]*c2 + Ay[2][2]*s2) + s2*(Ay[1][1]*c2 + Ay[1][2]*s2);
	R2AyR2[2][2] = c2*(Ay[2][2]*c2 - Ay[2][1]*s2) + s2*(Ay[1][2]*c2 - Ay[1][1]*s2);
	R2AyR2[2][3] = Ay[1][3]*s2 + Ay[2][3]*c2;
	R2AyR2[3][0] = Ay[3][0];
	R2AyR2[3][1] = Ay[3][2]*s2 + Ay[3][1]*c2;
	R2AyR2[3][2] = Ay[3][2]*c2 - Ay[3][1]*s2;
	R2AyR2[3][3] = Ay[3][3];

	//R2AyR3	
	R2AyR3[0][0] = Ay[0][0];
	R2AyR3[0][1] = Ay[0][2]*s3 + Ay[0][1]*c3;
	R2AyR3[0][2] = Ay[0][2]*c3 - Ay[0][1]*s3;
	R2AyR3[0][3] = Ay[0][3];
	R2AyR3[1][0] = Ay[1][0]*c2 - Ay[2][0]*s2;
	R2AyR3[1][1] = c2*(Ay[1][1]*c3 + Ay[1][2]*s3) - s2*(Ay[2][1]*c3 + Ay[2][2]*s3);
	R2AyR3[1][2] = c2*(Ay[1][2]*c3 - Ay[1][1]*s3) - s2*(Ay[2][2]*c3 - Ay[2][1]*s3);
	R2AyR3[1][3] = Ay[1][3]*c2 - Ay[2][3]*s2;
	R2AyR3[2][0] = Ay[1][0]*s2 + Ay[2][0]*c2;
	R2AyR3[2][1] = c2*(Ay[2][1]*c3 + Ay[2][2]*s3) + s2*(Ay[1][1]*c3 + Ay[1][2]*s3);
	R2AyR3[2][2] = c2*(Ay[2][2]*c3 - Ay[2][1]*s3) + s2*(Ay[1][2]*c3 - Ay[1][1]*s3);
	R2AyR3[2][3] = Ay[1][3]*s2 + Ay[2][3]*c2;
	R2AyR3[3][0] = Ay[3][0];
	R2AyR3[3][1] = Ay[3][2]*s3 + Ay[3][1]*c3;
	R2AyR3[3][2] = Ay[3][2]*c3 - Ay[3][1]*s3;
	R2AyR3[3][3] = Ay[3][3];

	//R3AyR1
	R3AyR1[0][0] = Ay[0][0];
	R3AyR1[0][1] = Ay[0][2]*s1 + Ay[0][1]*c1;
	R3AyR1[0][2] = Ay[0][2]*c1 - Ay[0][1]*s1;
	R3AyR1[0][3] = Ay[0][3];
	R3AyR1[1][0] = Ay[1][0]*c3 - Ay[2][0]*s3;
	R3AyR1[1][1] = c3*(Ay[1][1]*c1 + Ay[1][2]*s1) - s3*(Ay[2][1]*c1 + Ay[2][2]*s1);
	R3AyR1[1][2] = c3*(Ay[1][2]*c1 - Ay[1][1]*s1) - s3*(Ay[2][2]*c1 - Ay[2][1]*s1);
	R3AyR1[1][3] = Ay[1][3]*c3 - Ay[2][3]*s3;
	R3AyR1[2][0] = Ay[1][0]*s3 + Ay[2][0]*c3;
	R3AyR1[2][1] = c3*(Ay[2][1]*c1 + Ay[2][2]*s1) + s3*(Ay[1][1]*c1 + Ay[1][2]*s1);
	R3AyR1[2][2] = c3*(Ay[2][2]*c1 - Ay[2][1]*s1) + s3*(Ay[1][2]*c1 - Ay[1][1]*s1);
	R3AyR1[2][3] = Ay[1][3]*s3 + Ay[2][3]*c3;
	R3AyR1[3][0] = Ay[3][0];
	R3AyR1[3][1] = Ay[3][2]*s1 + Ay[3][1]*c1;
	R3AyR1[3][2] = Ay[3][2]*c1 - Ay[3][1]*s1;
	R3AyR1[3][3] = Ay[3][3];

	//R3AyR2	
	R3AyR2[0][0] = Ay[0][0];
	R3AyR2[0][1] = Ay[0][2]*s2 + Ay[0][1]*c2;
	R3AyR2[0][2] = Ay[0][2]*c2 - Ay[0][1]*s2;
	R3AyR2[0][3] = Ay[0][3];
	R3AyR2[1][0] = Ay[1][0]*c3 - Ay[2][0]*s3;
	R3AyR2[1][1] = c3*(Ay[1][1]*c2 + Ay[1][2]*s2) - s3*(Ay[2][1]*c2 + Ay[2][2]*s2);
	R3AyR2[1][2] = c3*(Ay[1][2]*c2 - Ay[1][1]*s2) - s3*(Ay[2][2]*c2 - Ay[2][1]*s2);
	R3AyR2[1][3] = Ay[1][3]*c3 - Ay[2][3]*s3;
	R3AyR2[2][0] = Ay[1][0]*s3 + Ay[2][0]*c3;
	R3AyR2[2][1] = c3*(Ay[2][1]*c2 + Ay[2][2]*s2) + s3*(Ay[1][1]*c2 + Ay[1][2]*s2);
	R3AyR2[2][2] = c3*(Ay[2][2]*c2 - Ay[2][1]*s2) + s3*(Ay[1][2]*c2 - Ay[1][1]*s2);
	R3AyR2[2][3] = Ay[1][3]*s3 + Ay[2][3]*c3;
	R3AyR2[3][0] = Ay[3][0];
	R3AyR2[3][1] = Ay[3][2]*s2 + Ay[3][1]*c2;
	R3AyR2[3][2] = Ay[3][2]*c2 - Ay[3][1]*s2;
	R3AyR2[3][3] = Ay[3][3];

	//R3AyR3	
	R3AyR3[0][0] = Ay[0][0];
	R3AyR3[0][1] = Ay[0][2]*s3 + Ay[0][1]*c3;
	R3AyR3[0][2] = Ay[0][2]*c3 - Ay[0][1]*s3;
	R3AyR3[0][3] = Ay[0][3];
	R3AyR3[1][0]= Ay[1][0]*c3 - Ay[2][0]*s3;
	R3AyR3[1][1]= c3*(Ay[1][1]*c3 + Ay[1][2]*s3) - s3*(Ay[2][1]*c3 + Ay[2][2]*s3);
	R3AyR3[1][2]= c3*(Ay[1][2]*c3 - Ay[1][1]*s3) - s3*(Ay[2][2]*c3 - Ay[2][1]*s3);
	R3AyR3[1][3]= Ay[1][3]*c3 - Ay[2][3]*s3;
	R3AyR3[2][0]= Ay[1][0]*s3 + Ay[2][0]*c3;
	R3AyR3[2][1]= c3*(Ay[2][1]*c3 + Ay[2][2]*s3) + s3*(Ay[1][1]*c3 + Ay[1][2]*s3);
	R3AyR3[2][2]= c3*(Ay[2][2]*c3 - Ay[2][1]*s3) + s3*(Ay[1][2]*c3 - Ay[1][1]*s3);
	R3AyR3[2][3]= Ay[1][3]*s3 + Ay[2][3]*c3;
	R3AyR3[3][0] = Ay[3][0];
	R3AyR3[3][1]= Ay[3][2]*s3 + Ay[3][1]*c3;
	R3AyR3[3][2]= Ay[3][2]*c3 - Ay[3][1]*s3;
	R3AyR3[3][3] = Ay[3][3];

	//R1Ay
	R1Ay[0][0] = Ay[0][0];
	R1Ay[0][1] = Ay[0][1];
	R1Ay[0][2] = Ay[0][2];
	R1Ay[0][3] = Ay[0][3];
	R1Ay[1][0] = Ay[1][0]*c1 - Ay[2][0]*s1;
	R1Ay[1][1] = Ay[1][1]*c1 - Ay[2][1]*s1;
	R1Ay[1][2] = Ay[1][2]*c1 - Ay[2][2]*s1;
	R1Ay[1][3] = Ay[1][3]*c1 - Ay[2][3]*s1;
	R1Ay[2][0] = Ay[1][0]*s1 + Ay[2][0]*c1;
	R1Ay[2][1] = Ay[1][1]*s1 + Ay[2][1]*c1;
	R1Ay[2][2] = Ay[1][2]*s1 + Ay[2][2]*c1;
	R1Ay[2][3] = Ay[1][3]*s1 + Ay[2][3]*c1;
	R1Ay[3][0] = Ay[3][0];
	R1Ay[3][1] = Ay[3][1];
	R1Ay[3][2] = Ay[3][2];
	R1Ay[3][3] = Ay[3][3];
	
	//R2Ay
	R2Ay[0][0] = Ay[0][0];
	R2Ay[0][1] = Ay[0][1];
	R2Ay[0][2] = Ay[0][2];
	R2Ay[0][3] = Ay[0][3];
	R2Ay[1][0] = Ay[1][0]*c2 - Ay[2][0]*s2;
	R2Ay[1][1] = Ay[1][1]*c2 - Ay[2][1]*s2;
	R2Ay[1][2] = Ay[1][2]*c2 - Ay[2][2]*s2;
	R2Ay[1][3] = Ay[1][3]*c2 - Ay[2][3]*s2;
	R2Ay[2][0] = Ay[1][0]*s2 + Ay[2][0]*c2;
	R2Ay[2][1] = Ay[1][1]*s2 + Ay[2][1]*c2;
	R2Ay[2][2] = Ay[1][2]*s2 + Ay[2][2]*c2;
	R2Ay[2][3] = Ay[1][3]*s2 + Ay[2][3]*c2;
	R2Ay[3][0] = Ay[3][0];
	R2Ay[3][1] = Ay[3][1];
	R2Ay[3][2] = Ay[3][2];
	R2Ay[3][3] = Ay[3][3];
	
	//R3Ay
	R3Ay[0][0] = Ay[0][0];
	R3Ay[0][1] = Ay[0][1];
	R3Ay[0][2] = Ay[0][2];
	R3Ay[0][3] = Ay[0][3];
	R3Ay[1][0] = Ay[1][0]*c3 - Ay[2][0]*s3;
	R3Ay[1][1] = Ay[1][1]*c3 - Ay[2][1]*s3;
	R3Ay[1][2] = Ay[1][2]*c3 - Ay[2][2]*s3;
	R3Ay[1][3] = Ay[1][3]*c3 - Ay[2][3]*s3;
	R3Ay[2][0] = Ay[1][0]*s3 + Ay[2][0]*c3;
	R3Ay[2][1] = Ay[1][1]*s3 + Ay[2][1]*c3;
	R3Ay[2][2] = Ay[1][2]*s3 + Ay[2][2]*c3;
	R3Ay[2][3] = Ay[1][3]*s3 + Ay[2][3]*c3;
	R3Ay[3][0] = Ay[3][0];
	R3Ay[3][1] = Ay[3][1];
	R3Ay[3][2] = Ay[3][2];
	R3Ay[3][3] = Ay[3][3];

	//AyR1
	AyR1[0][0] = Ay[0][0];
	AyR1[0][1] = Ay[0][2]*s1 + Ay[0][1]*c1;
	AyR1[0][2] = Ay[0][2]*c1 - Ay[0][1]*s1;
	AyR1[0][3] = Ay[0][3];	
	AyR1[1][0] = Ay[1][0];	
	AyR1[1][1] = Ay[1][2]*s1 + Ay[1][1]*c1;
	AyR1[1][2] = Ay[1][2]*c1 - Ay[1][1]*s1;
	AyR1[1][3] = Ay[1][3];	
	AyR1[2][0] = Ay[2][0];	
	AyR1[2][1] = Ay[2][2]*s1 + Ay[2][1]*c1;
	AyR1[2][2] = Ay[2][2]*c1 - Ay[2][1]*s1;
	AyR1[2][3] = Ay[2][3];	
	AyR1[3][0] = Ay[3][0];	
	AyR1[3][1] = Ay[3][2]*s1 + Ay[3][1]*c1;
	AyR1[3][2] = Ay[3][2]*c1 - Ay[3][1]*s1;
	AyR1[3][3] = Ay[3][3];	
	
	//AyR2
	AyR2[0][0] = Ay[0][0];
	AyR2[0][1] = Ay[0][2]*s2 + Ay[0][1]*c2;
	AyR2[0][2] = Ay[0][2]*c2 - Ay[0][1]*s2;
	AyR2[0][3] = Ay[0][3];	
	AyR2[1][0] = Ay[1][0];	
	AyR2[1][1] = Ay[1][2]*s2 + Ay[1][1]*c2;
	AyR2[1][2] = Ay[1][2]*c2 - Ay[1][1]*s2;
	AyR2[1][3] = Ay[1][3];	
	AyR2[2][0] = Ay[2][0];	
	AyR2[2][1] = Ay[2][2]*s2 + Ay[2][1]*c2;
	AyR2[2][2] = Ay[2][2]*c2 - Ay[2][1]*s2;
	AyR2[2][3] = Ay[2][3];	
	AyR2[3][0] = Ay[3][0];	
	AyR2[3][1] = Ay[3][2]*s2 + Ay[3][1]*c2;
	AyR2[3][2] = Ay[3][2]*c2 - Ay[3][1]*s2;
	AyR2[3][3] = Ay[3][3];	

	//AyR3
	AyR3[0][0] = Ay[0][0];
	AyR3[0][1] = Ay[0][2]*s3 + Ay[0][1]*c3;
	AyR3[0][2] = Ay[0][2]*c3 - Ay[0][1]*s3;
	AyR3[0][3] = Ay[0][3];	
	AyR3[1][0] = Ay[1][0];	
	AyR3[1][1] = Ay[1][2]*s3 + Ay[1][1]*c3;
	AyR3[1][2] = Ay[1][2]*c3 - Ay[1][1]*s3;
	AyR3[1][3] = Ay[1][3];	
	AyR3[2][0] = Ay[2][0];	
	AyR3[2][1] = Ay[2][2]*s3 + Ay[2][1]*c3;
	AyR3[2][2] = Ay[2][2]*c3 - Ay[2][1]*s3;
	AyR3[2][3] = Ay[2][3];	
	AyR3[3][0] = Ay[3][0];	
	AyR3[3][1] = Ay[3][2]*s3 + Ay[3][1]*c3;
	AyR3[3][2] = Ay[3][2]*c3 - Ay[3][1]*s3;
	AyR3[3][3] = Ay[3][3];	

	//R1R1
	R1R1[0][0] = 1;
	R1R1[0][1] = 0;
	R1R1[0][2] = 0;
	R1R1[0][3] = 0;
	R1R1[1][0] = 0;
	R1R1[1][1] =  c1*c1 - s1*s1;
	R1R1[1][2] = -c1*s1 - c1*s1;
	R1R1[1][3] = 0;
	R1R1[2][0] = 0;
	R1R1[2][1] =  c1*s1 + c1*s1;
	R1R1[2][2] =  c1*c1 - s1*s1;
	R1R1[2][3] = 0;
	R1R1[3][0] = 0;
	R1R1[3][1] = 0;
	R1R1[3][2] = 0;
	R1R1[3][3] = 1;

	//R1R2
	R1R2[0][0] = 1;
	R1R2[0][1] = 0;
	R1R2[0][2] = 0;
	R1R2[0][3] = 0;
	R1R2[1][0] = 0;
	R1R2[1][1] =  c2*c1 - s2*s1;
	R1R2[1][2] = -c2*s1 - c1*s2;
	R1R2[1][3] = 0;
	R1R2[2][0] = 0;
	R1R2[2][1] =  c2*s1 + c1*s2;
	R1R2[2][2] =  c2*c1 - s2*s1;
	R1R2[2][3] = 0;
	R1R2[3][0] = 0;
	R1R2[3][1] = 0;
	R1R2[3][2] = 0;
	R1R2[3][3] = 1;

	//R1R3
	R1R3[0][0] = 1;
	R1R3[0][1] = 0;
	R1R3[0][2] = 0;
	R1R3[0][3] = 0;
	R1R3[1][0] = 0;
	R1R3[1][1] =  c3*c1 - s3*s1;
	R1R3[1][2] = -c3*s1 - c1*s3;
	R1R3[1][3] = 0;
	R1R3[2][0] = 0;
	R1R3[2][1] =  c3*s1 + c1*s3;
	R1R3[2][2] =  c3*c1 - s3*s1;
	R1R3[2][3] = 0;
	R1R3[3][0] = 0;
	R1R3[3][1] = 0;
	R1R3[3][2] = 0;
	R1R3[3][3] = 1;
	
	//R2R2
	R2R2[0][0] = 1;
	R2R2[0][1] = 0;
	R2R2[0][2] = 0;
	R2R2[0][3] = 0;
	R2R2[1][0] = 0;
	R2R2[1][1] =  c2*c2 - s2*s2;
	R2R2[1][2] = -c2*s2 - c2*s2;
	R2R2[1][3] = 0;
	R2R2[2][0] = 0;
	R2R2[2][1] =  c2*s2 + c2*s2;
	R2R2[2][2] =  c2*c2 - s2*s2;
	R2R2[2][3] = 0;
	R2R2[3][0] = 0;
	R2R2[3][1] = 0;
	R2R2[3][2] = 0;
	R2R2[3][3] = 1;

	//R2R3
	R2R3[0][0] = 1;
	R2R3[0][1] = 0;
	R2R3[0][2] = 0;
	R2R3[0][3] = 0;
	R2R3[1][0] = 0;
	R2R3[1][1] =  c3*c2 - s3*s2;
	R2R3[1][2] = -c3*s2 - c1*s3;
	R2R3[1][3] = 0;
	R2R3[2][0] = 0;
	R2R3[2][1] =  c3*s2 + c2*s3;
	R2R3[2][2] =  c3*c2 - s3*s2;
	R2R3[2][3] = 0;
	R2R3[3][0] = 0;
	R2R3[3][1] = 0;
	R2R3[3][2] = 0;
	R2R3[3][3] = 1;

	//R3R3
	R3R3[0][0] = 1;
	R3R3[0][1] = 0;
	R3R3[0][2] = 0;
	R3R3[0][3] = 0;
	R3R3[1][0] = 0;
	R3R3[1][1] =  c3*c3 - s3*s3;
	R3R3[1][2] = -c3*s3 - c3*s3;
	R3R3[1][3] = 0;
	R3R3[2][0] = 0;
	R3R3[2][1] =  c3*s3 + c3*s3;
	R3R3[2][2] =  c3*c3 - s3*s3;
	R3R3[2][3] = 0;
	R3R3[3][0] = 0;
	R3R3[3][1] = 0;
	R3R3[3][2] = 0;
	R3R3[3][3] = 1;

	//R1
	R1[0][0] = 1;
	R1[0][1] = 0;
	R1[0][2] = 0;
	R1[0][3] = 0;
	R1[1][0] = 0;
	R1[1][1] = c1;
	R1[1][2] = -s1;
	R1[1][3] = 0;
	R1[2][0] = 0;
	R1[2][1] = s1;
	R1[2][2] = c1;
	R1[2][3] = 0;
	R1[3][0] = 0;
	R1[3][1] = 0;
	R1[3][2] = 0;
	R1[3][3] = 1;
	
	//R2
	R2[0][0] = 1;
	R2[0][1] = 0;
	R2[0][2] = 0;
	R2[0][3] = 0;
	R2[1][0] = 0;
	R2[1][1] = c2;
	R2[1][2] = -s2;
	R2[1][3] = 0;
	R2[2][0] = 0;
	R2[2][1] = s2;
	R2[2][2] = c2;
	R2[2][3] = 0;
	R2[3][0] = 0;
	R2[3][1] = 0;
	R2[3][2] = 0;
	R2[3][3] = 1;
	
	//R3
	R3[0][0] = 1;
	R3[0][1] = 0;
	R3[0][2] = 0;
	R3[0][3] = 0;
	R3[1][0] = 0;
	R3[1][1] = c3;
	R3[1][2] = -s3;
	R3[1][3] = 0;
	R3[2][0] = 0;
	R3[2][1] = s3;
	R3[2][2] = c3;
	R3[2][3] = 0;
	R3[3][0] = 0;
	R3[3][1] = 0;
	R3[3][2] = 0;
	R3[3][3] = 1;


	//Bloco K11:
	double k12, k13, k21, k22, k23, k24, k31, k32, k33, k34, k42, k43;
	
 	k12 = one_sixth*(y23*R1AxR1[0][1] + x32*R1AyR1[0][1]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[0][1];
 	k13 = one_sixth*(y23*R1AxR1[0][2] + x32*R1AyR1[0][2]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[0][2];
 	k21 = one_sixth*(y23*R1AxR1[1][0] + x32*R1AyR1[1][0]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[1][0];
 	k22 = one_sixth*(y23*R1AxR1[1][1] + x32*R1AyR1[1][1]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[1][1];
 	k23 = one_sixth*(y23*R1AxR1[1][2] + x32*R1AyR1[1][2]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[1][2];
 	k24 = one_sixth*(y23*R1AxR1[1][3] + x32*R1AyR1[1][3]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[1][3];
 	k31 = one_sixth*(y23*R1AxR1[2][0] + x32*R1AyR1[2][0]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[2][0];
 	k32 = one_sixth*(y23*R1AxR1[2][1] + x32*R1AyR1[2][1]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[2][1];
 	k33 = one_sixth*(y23*R1AxR1[2][2] + x32*R1AyR1[2][2]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[2][2];
 	k34 = one_sixth*(y23*R1AxR1[2][3] + x32*R1AyR1[2][3]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[2][3];
 	k42 = one_sixth*(y23*R1AxR1[3][1] + x32*R1AyR1[3][1]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[3][1];
 	k43 = one_sixth*(y23*R1AxR1[3][2] + x32*R1AyR1[3][2]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[3][2];
 	
	//Bloco K12:
	double k16, k17, k25, k26, k27, k28, k35, k36, k37, k38, k46, k47; 
	k16 = one_sixth*(y31*R1AxR2[0][1] + x13*R1AyR2[0][1]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[0][1];
 	k17 = one_sixth*(y31*R1AxR2[0][2] + x13*R1AyR2[0][2]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[0][2];
 	k25 = one_sixth*(y31*R1AxR2[1][0] + x13*R1AyR2[1][0]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[1][0];
 	k26 = one_sixth*(y31*R1AxR2[1][1] + x13*R1AyR2[1][1]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[1][1];
 	k27 = one_sixth*(y31*R1AxR2[1][2] + x13*R1AyR2[1][2]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[1][2];
 	k28 = one_sixth*(y31*R1AxR2[1][3] + x13*R1AyR2[1][3]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[1][3];
 	k35 = one_sixth*(y31*R1AxR2[2][0] + x13*R1AyR2[2][0]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[2][0];
 	k36 = one_sixth*(y31*R1AxR2[2][1] + x13*R1AyR2[2][1]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[2][1];
 	k37 = one_sixth*(y31*R1AxR2[2][2] + x13*R1AyR2[2][2]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[2][2];
 	k38 = one_sixth*(y31*R1AxR2[2][3] + x13*R1AyR2[2][3]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[2][3];
 	k46 = one_sixth*(y31*R1AxR2[3][1] + x13*R1AyR2[3][1]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[3][1];
 	k47 = one_sixth*(y31*R1AxR2[3][2] + x13*R1AyR2[3][2]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[3][2];

	//Bloco K13
	double k110, k1_11, k29, k210, k211, k212, k39, k310, k311, k312, k410, k411;
	k110 = one_sixth*(y12*R1AxR3[0][1] + x21*R1AyR3[0][1]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[0][1];
 	k1_11 = one_sixth*(y12*R1AxR3[0][2] + x21*R1AyR3[0][2]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[0][2];
 	k29  = one_sixth*(y12*R1AxR3[1][0] + x21*R1AyR3[1][0]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[1][0];
 	k210 = one_sixth*(y12*R1AxR3[1][1] + x21*R1AyR3[1][1]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[1][1];
 	k211 = one_sixth*(y12*R1AxR3[1][2] + x21*R1AyR3[1][2]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[1][2];
 	k212 = one_sixth*(y12*R1AxR3[1][3] + x21*R1AyR3[1][3]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[1][3];
 	k39  = one_sixth*(y12*R1AxR3[2][0] + x21*R1AyR3[2][0]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[2][0];
 	k310 = one_sixth*(y12*R1AxR3[2][1] + x21*R1AyR3[2][1]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[2][1];
 	k311 = one_sixth*(y12*R1AxR3[2][2] + x21*R1AyR3[2][2]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[2][2];
 	k312 = one_sixth*(y12*R1AxR3[2][3] + x21*R1AyR3[2][3]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[2][3];
 	k410 = one_sixth*(y12*R1AxR3[3][1] + x21*R1AyR3[3][1]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[3][1];
 	k411 = one_sixth*(y12*R1AxR3[3][2] + x21*R1AyR3[3][2]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[3][2];

 	//Bloco K14
	double k213, k214, k215, k216, k313, k314, k315, k316;

	k213 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[1][0] + (x32 + 2*x13 + 2*x21)*R1Ay[1][0]);		
	k214 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[1][1] + (x32 + 2*x13 + 2*x21)*R1Ay[1][1]);		
	k215 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[1][2] + (x32 + 2*x13 + 2*x21)*R1Ay[1][2]);		
	k216 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[1][3] + (x32 + 2*x13 + 2*x21)*R1Ay[1][3]);		
	k313 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[2][0] + (x32 + 2*x13 + 2*x21)*R1Ay[2][0]);		
	k314 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[2][1] + (x32 + 2*x13 + 2*x21)*R1Ay[2][1]);		
	k315 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[2][2] + (x32 + 2*x13 + 2*x21)*R1Ay[2][2]);		
	k316 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[2][3] + (x32 + 2*x13 + 2*x21)*R1Ay[2][3]);		

	//Bloco K21
	double k52, k53, k61, k62, k63, k64, k71, k72, k73, k74, k82, k83;

	k52 = one_sixth*(y23*R2AxR1[0][1] + x32*R2AyR1[0][1]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[0][1];
 	k53 = one_sixth*(y23*R2AxR1[0][2] + x32*R2AyR1[0][2]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[0][2];
 	k61 = one_sixth*(y23*R2AxR1[1][0] + x32*R2AyR1[1][0]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[1][0];
 	k62 = one_sixth*(y23*R2AxR1[1][1] + x32*R2AyR1[1][1]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[1][1];
 	k63 = one_sixth*(y23*R2AxR1[1][2] + x32*R2AyR1[1][2]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[1][2];
 	k64 = one_sixth*(y23*R2AxR1[1][3] + x32*R2AyR1[1][3]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[1][3];
 	k71 = one_sixth*(y23*R2AxR1[2][0] + x32*R2AyR1[2][0]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[2][0];
 	k72 = one_sixth*(y23*R2AxR1[2][1] + x32*R2AyR1[2][1]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[2][1];
 	k73 = one_sixth*(y23*R2AxR1[2][2] + x32*R2AyR1[2][2]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[2][2];
 	k74 = one_sixth*(y23*R2AxR1[2][3] + x32*R2AyR1[2][3]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[2][3];
 	k82 = one_sixth*(y23*R2AxR1[3][1] + x32*R2AyR1[3][1]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[3][1];
 	k83 = one_sixth*(y23*R2AxR1[3][2] + x32*R2AyR1[3][2]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[3][2];

	
	//Bloco K22
	double k56, k57, k65, k66, k67, k68, k75, k76, k77, k78, k86, k87;
	
	k56 = one_sixth*(y31*R2AxR2[0][1] + x13*R2AyR2[0][1]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[0][1];
 	k57 = one_sixth*(y31*R2AxR2[0][2] + x13*R2AyR2[0][2]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[0][2];
 	k65 = one_sixth*(y31*R2AxR2[1][0] + x13*R2AyR2[1][0]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[1][0];
 	k66 = one_sixth*(y31*R2AxR2[1][1] + x13*R2AyR2[1][1]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[1][1];
 	k67 = one_sixth*(y31*R2AxR2[1][2] + x13*R2AyR2[1][2]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[1][2];
 	k68 = one_sixth*(y31*R2AxR2[1][3] + x13*R2AyR2[1][3]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[1][3];
 	k75 = one_sixth*(y31*R2AxR2[2][0] + x13*R2AyR2[2][0]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[2][0];
 	k76 = one_sixth*(y31*R2AxR2[2][1] + x13*R2AyR2[2][1]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[2][1];
 	k77 = one_sixth*(y31*R2AxR2[2][2] + x13*R2AyR2[2][2]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[2][2];
 	k78 = one_sixth*(y31*R2AxR2[2][3] + x13*R2AyR2[2][3]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[2][3];
 	k86 = one_sixth*(y31*R2AxR2[3][1] + x13*R2AyR2[3][1]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[3][1];
 	k87 = one_sixth*(y31*R2AxR2[3][2] + x13*R2AyR2[3][2]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[3][2];

	//Bloco K23:
	double k510, k511, k69, k610, k611, k612, k79, k710, k711, k712, k810, k811;

	k510 = one_sixth*(y12*R2AxR3[0][1] + x21*R2AyR3[0][1]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[0][1];
 	k511 = one_sixth*(y12*R2AxR3[0][2] + x21*R2AyR3[0][2]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[0][2];
 	k69  = one_sixth*(y12*R2AxR3[1][0] + x21*R2AyR3[1][0]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[1][0];
 	k610 = one_sixth*(y12*R2AxR3[1][1] + x21*R2AyR3[1][1]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[1][1];
 	k611 = one_sixth*(y12*R2AxR3[1][2] + x21*R2AyR3[1][2]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[1][2];
 	k612 = one_sixth*(y12*R2AxR3[1][3] + x21*R2AyR3[1][3]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[1][3];
 	k79  = one_sixth*(y12*R2AxR3[2][0] + x21*R2AyR3[2][0]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[2][0];
 	k710 = one_sixth*(y12*R2AxR3[2][1] + x21*R2AyR3[2][1]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[2][1];
 	k711 = one_sixth*(y12*R2AxR3[2][2] + x21*R2AyR3[2][2]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[2][2];
 	k712 = one_sixth*(y12*R2AxR3[2][3] + x21*R2AyR3[2][3]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[2][3];
 	k810 = one_sixth*(y12*R2AxR3[3][1] + x21*R2AyR3[3][1]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[3][1];
 	k811 = one_sixth*(y12*R2AxR3[3][2] + x21*R2AyR3[3][2]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[3][2];

	//Bloco K24
	double k613, k614, k615, k616, k713, k714, k715, k716;

	k613 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[1][0] + (2*x32 + x13 + 2*x21)*R2Ay[1][0]);		
	k614 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[1][1] + (2*x32 + x13 + 2*x21)*R2Ay[1][1]);		
	k615 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[1][2] + (2*x32 + x13 + 2*x21)*R2Ay[1][2]);		
	k616 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[1][3] + (2*x32 + x13 + 2*x21)*R2Ay[1][3]);		
	k713 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[2][0] + (2*x32 + x13 + 2*x21)*R2Ay[2][0]);		
	k714 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[2][1] + (2*x32 + x13 + 2*x21)*R2Ay[2][1]);		
	k715 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[2][2] + (2*x32 + x13 + 2*x21)*R2Ay[2][2]);		
	k716 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[2][3] + (2*x32 + x13 + 2*x21)*R2Ay[2][3]);		

	//Bloco K31
	double k92, k93, k101, k102, k103, k104, k11_1, k11_2, k113, k114, k122, k123;
	
	k92   = one_sixth*(y23*R3AxR1[0][1] + x32*R3AyR1[0][1]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[0][1];
 	k93   = one_sixth*(y23*R3AxR1[0][2] + x32*R3AyR1[0][2]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[0][2];
 	k101  = one_sixth*(y23*R3AxR1[1][0] + x32*R3AyR1[1][0]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[1][0];
 	k102  = one_sixth*(y23*R3AxR1[1][1] + x32*R3AyR1[1][1]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[1][1];
 	k103  = one_sixth*(y23*R3AxR1[1][2] + x32*R3AyR1[1][2]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[1][2];
 	k104  = one_sixth*(y23*R3AxR1[1][3] + x32*R3AyR1[1][3]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[1][3];
 	k11_1 = one_sixth*(y23*R3AxR1[2][0] + x32*R3AyR1[2][0]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[2][0];
 	k11_2 = one_sixth*(y23*R3AxR1[2][1] + x32*R3AyR1[2][1]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[2][1];
 	k113  = one_sixth*(y23*R3AxR1[2][2] + x32*R3AyR1[2][2]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[2][2];
 	k114  = one_sixth*(y23*R3AxR1[2][3] + x32*R3AyR1[2][3]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[2][3];
 	k122  = one_sixth*(y23*R3AxR1[3][1] + x32*R3AyR1[3][1]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[3][1];
 	k123  = one_sixth*(y23*R3AxR1[3][2] + x32*R3AyR1[3][2]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[3][2];

	
	//Bloco K32
	double k96, k97, k105, k106, k107, k108, k115, k116, k117, k118, k126, k127;

	k96  = one_sixth*(y31*R3AxR2[0][1] + x13*R3AyR2[0][1]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[0][1];
 	k97  = one_sixth*(y31*R3AxR2[0][2] + x13*R3AyR2[0][2]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[0][2];
 	k105 = one_sixth*(y31*R3AxR2[1][0] + x13*R3AyR2[1][0]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[1][0];
 	k106 = one_sixth*(y31*R3AxR2[1][1] + x13*R3AyR2[1][1]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[1][1];
 	k107 = one_sixth*(y31*R3AxR2[1][2] + x13*R3AyR2[1][2]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[1][2];
 	k108 = one_sixth*(y31*R3AxR2[1][3] + x13*R3AyR2[1][3]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[1][3];
 	k115 = one_sixth*(y31*R3AxR2[2][0] + x13*R3AyR2[2][0]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[2][0];
 	k116 = one_sixth*(y31*R3AxR2[2][1] + x13*R3AyR2[2][1]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[2][1];
 	k117 = one_sixth*(y31*R3AxR2[2][2] + x13*R3AyR2[2][2]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[2][2];
 	k118 = one_sixth*(y31*R3AxR2[2][3] + x13*R3AyR2[2][3]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[2][3];
 	k126 = one_sixth*(y31*R3AxR2[3][1] + x13*R3AyR2[3][1]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[3][1];
 	k127 = one_sixth*(y31*R3AxR2[3][2] + x13*R3AyR2[3][2]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[3][2];

	//Bloco K33
	double k910, k911, k109, k1010, k1011, k1012, k119, k1110, k1111, k1112, k1210, k1211;
	
	k910  = one_sixth*(y12*R3AxR3[0][1] + x21*R3AyR3[0][1]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[0][1];
 	k911  = one_sixth*(y12*R3AxR3[0][2] + x21*R3AyR3[0][2]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[0][2];
 	k109  = one_sixth*(y12*R3AxR3[1][0] + x21*R3AyR3[1][0]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[1][0];
 	k1010 = one_sixth*(y12*R3AxR3[1][1] + x21*R3AyR3[1][1]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[1][1];
 	k1011 = one_sixth*(y12*R3AxR3[1][2] + x21*R3AyR3[1][2]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[1][2];
 	k1012 = one_sixth*(y12*R3AxR3[1][3] + x21*R3AyR3[1][3]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[1][3];
 	k119  = one_sixth*(y12*R3AxR3[2][0] + x21*R3AyR3[2][0]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[2][0];
 	k1110 = one_sixth*(y12*R3AxR3[2][1] + x21*R3AyR3[2][1]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[2][1];
 	k1111 = one_sixth*(y12*R3AxR3[2][2] + x21*R3AyR3[2][2]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[2][2];
 	k1112 = one_sixth*(y12*R3AxR3[2][3] + x21*R3AyR3[2][3]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[2][3];
 	k1210 = one_sixth*(y12*R3AxR3[3][1] + x21*R3AyR3[3][1]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[3][1];
 	k1211 = one_sixth*(y12*R3AxR3[3][2] + x21*R3AyR3[3][2]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[3][2];


	//Bloco K34
	double k1013, k1014, k1015, k1016, k1113, k1114, k1115, k1116;

	k1013 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[1][0] + (2*x32 + 2*x13 + x21)*R3Ay[1][0]);		
	k1014 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[1][1] + (2*x32 + 2*x13 + x21)*R3Ay[1][1]);		
	k1015 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[1][2] + (2*x32 + 2*x13 + x21)*R3Ay[1][2]);		
	k1016 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[1][3] + (2*x32 + 2*x13 + x21)*R3Ay[1][3]);		
	k1113 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[2][0] + (2*x32 + 2*x13 + x21)*R3Ay[2][0]);		
	k1114 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[2][1] + (2*x32 + 2*x13 + x21)*R3Ay[2][1]);		
	k1115 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[2][2] + (2*x32 + 2*x13 + x21)*R3Ay[2][2]);		
	k1116 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[2][3] + (2*x32 + 2*x13 + x21)*R3Ay[2][3]);		

	//Bloco K41:
	double k132, k133, k142, k143, k152, k153, k162, k163;

	k132 = nine_over_forty*(y23*AxR1[0][1] + x32*AyR1[0][1]);
	k133 = nine_over_forty*(y23*AxR1[0][2] + x32*AyR1[0][2]);
	k142 = nine_over_forty*(y23*AxR1[1][1] + x32*AyR1[1][1]);
	k143 = nine_over_forty*(y23*AxR1[1][2] + x32*AyR1[1][2]);
	k152 = nine_over_forty*(y23*AxR1[2][1] + x32*AyR1[2][1]);
	k153 = nine_over_forty*(y23*AxR1[2][2] + x32*AyR1[2][2]);
	k162 = nine_over_forty*(y23*AxR1[3][1] + x32*AyR1[3][1]);
	k163 = nine_over_forty*(y23*AxR1[3][2] + x32*AyR1[3][2]);

	//Bloco K42:
	double k136, k137, k146, k147, k156, k157, k166, k167;

	k136 = nine_over_forty*(y31*AxR2[0][1] + x13*AyR2[0][1]);
	k137 = nine_over_forty*(y31*AxR2[0][2] + x13*AyR2[0][2]);
	k146 = nine_over_forty*(y31*AxR2[1][1] + x13*AyR2[1][1]);
	k147 = nine_over_forty*(y31*AxR2[1][2] + x13*AyR2[1][2]);
	k156 = nine_over_forty*(y31*AxR2[2][1] + x13*AyR2[2][1]);
	k157 = nine_over_forty*(y31*AxR2[2][2] + x13*AyR2[2][2]);
	k166 = nine_over_forty*(y31*AxR2[3][1] + x13*AyR2[3][1]);
	k167 = nine_over_forty*(y31*AxR2[3][2] + x13*AyR2[3][2]);
	
	//Bloco K43:
	double k1310, k1311, k1410, k1411, k1510, k1511, k1610, k1611;

	k1310 = nine_over_forty*(y12*AxR3[0][1] + x21*AyR3[0][1]);
	k1311 = nine_over_forty*(y12*AxR3[0][2] + x21*AyR3[0][2]);
	k1410 = nine_over_forty*(y12*AxR3[1][1] + x21*AyR3[1][1]);
	k1411 = nine_over_forty*(y12*AxR3[1][2] + x21*AyR3[1][2]);
	k1510 = nine_over_forty*(y12*AxR3[2][1] + x21*AyR3[2][1]);
	k1511 = nine_over_forty*(y12*AxR3[2][2] + x21*AyR3[2][2]);
	k1610 = nine_over_forty*(y12*AxR3[3][1] + x21*AyR3[3][1]);
	k1611 = nine_over_forty*(y12*AxR3[3][2] + x21*AyR3[3][2]);
	
	//Block M11;		
	double m12, m13, m21, m22, m23, m24, m31, m32, m33, m34, m42, m43;

	m12 = Area_over_six*R1R1[0][1];
	m13 = Area_over_six*R1R1[0][2];
	m21 = Area_over_six*R1R1[1][0];
	m22 = Area_over_six*R1R1[1][1];
	m23 = Area_over_six*R1R1[1][2];
	m24 = Area_over_six*R1R1[1][3];
	m31 = Area_over_six*R1R1[2][0];
	m32 = Area_over_six*R1R1[2][1];
	m33 = Area_over_six*R1R1[2][2];
	m34 = Area_over_six*R1R1[2][3];
	m42 = Area_over_six*R1R1[3][1];
	m43 = Area_over_six*R1R1[3][2];

	//Bloco M12
	double m16, m17, m25, m26, m27, m28, m35, m36, m37, m38, m46, m47; 
	
	m16 = Area_over_twelve*R1R2[0][1];
	m17 = Area_over_twelve*R1R2[0][2];
	m25 = Area_over_twelve*R1R2[1][0];
	m26 = Area_over_twelve*R1R2[1][1];
	m27 = Area_over_twelve*R1R2[1][2];
	m28 = Area_over_twelve*R1R2[1][3];
	m35 = Area_over_twelve*R1R2[2][0];
	m36 = Area_over_twelve*R1R2[2][1];
	m37 = Area_over_twelve*R1R2[2][2];
	m38 = Area_over_twelve*R1R2[2][3];
	m46 = Area_over_twelve*R1R2[3][1];
	m47 = Area_over_twelve*R1R2[3][2];

	//Bloco M13
	double m110, m111, m29, m210, m211, m212, m39, m310, m311, m312, m410, m411;

	m110 = Area_over_twelve*R1R3[0][1];
	m111 = Area_over_twelve*R1R3[0][2];
	m29  = Area_over_twelve*R1R3[1][0];
	m210 = Area_over_twelve*R1R3[1][1];
	m211 = Area_over_twelve*R1R3[1][2];
	m212 = Area_over_twelve*R1R3[1][3];
	m39  = Area_over_twelve*R1R3[2][0];
	m310 = Area_over_twelve*R1R3[2][1];
	m311 = Area_over_twelve*R1R3[2][2];
	m312 = Area_over_twelve*R1R3[2][3];
	m410 = Area_over_twelve*R1R3[3][1];
	m411 = Area_over_twelve*R1R3[3][2];

	//Bloco M14
	double m213, m214, m215, m216, m313, m314, m315, m316;

	m213 = Area_3_over_twenty*R1[1][0];	
	m214 = Area_3_over_twenty*R1[1][1];	
	m215 = Area_3_over_twenty*R1[1][2];	
	m216 = Area_3_over_twenty*R1[1][3];	
	m313 = Area_3_over_twenty*R1[2][0];	
	m314 = Area_3_over_twenty*R1[2][1];	
	m315 = Area_3_over_twenty*R1[2][2];	
	m316 = Area_3_over_twenty*R1[2][3];	
	
	//Bloco M22
	double m56, m57, m65, m66, m67, m68, m75, m76, m77, m78, m86, m87;

	m56 = Area_over_six*R2R2[0][1];
	m57 = Area_over_six*R2R2[0][2];
	m65 = Area_over_six*R2R2[1][0];
	m66 = Area_over_six*R2R2[1][1];
	m67 = Area_over_six*R2R2[1][2];
	m68 = Area_over_six*R2R2[1][3];
	m75 = Area_over_six*R2R2[2][0];
	m76 = Area_over_six*R2R2[2][1];
	m77 = Area_over_six*R2R2[2][2];
	m78 = Area_over_six*R2R2[2][3];
	m86 = Area_over_six*R2R2[3][1];
	m87 = Area_over_six*R2R2[3][2];


	//Bloco M23
	double m510, m511, m69, m610, m611, m612, m79, m710, m711, m712, m810, m811;
	
	m510 = Area_over_twelve*R2R3[0][1];
	m511 = Area_over_twelve*R2R3[0][2];
	m69  = Area_over_twelve*R2R3[1][0];
	m610 = Area_over_twelve*R2R3[1][1];
	m611 = Area_over_twelve*R2R3[1][2];
	m612 = Area_over_twelve*R2R3[1][3];
	m79  = Area_over_twelve*R2R3[2][0];
	m710 = Area_over_twelve*R2R3[2][1];
	m711 = Area_over_twelve*R2R3[2][2];
	m712 = Area_over_twelve*R2R3[2][3];
	m810 = Area_over_twelve*R2R3[3][1];
	m811 = Area_over_twelve*R2R3[3][2];

	//Bloco M24
	double m613, m614, m615, m616, m713, m714, m715, m716;
	
	m613 = Area_3_over_twenty*R2[1][0];	
	m614 = Area_3_over_twenty*R2[1][1];	
	m615 = Area_3_over_twenty*R2[1][2];	
	m616 = Area_3_over_twenty*R2[1][3];	
	m713 = Area_3_over_twenty*R2[2][0];	
	m714 = Area_3_over_twenty*R2[2][1];	
	m715 = Area_3_over_twenty*R2[2][2];	
	m716 = Area_3_over_twenty*R2[2][3];	

	//Bloco M33
	double m910, m911, m109, m1010, m1011, m1012, m119, m1110, m1111, m1112, m1210, m1211;
	
	m910  = Area_over_six*R3R3[0][1];
	m911  = Area_over_six*R3R3[0][2];
	m109  = Area_over_six*R3R3[1][0];
	m1010 = Area_over_six*R3R3[1][1];
	m1011 = Area_over_six*R3R3[1][2];
	m1012 = Area_over_six*R3R3[1][3];
	m119  = Area_over_six*R3R3[2][0];
	m1110 = Area_over_six*R3R3[2][1];
	m1111 = Area_over_six*R3R3[2][2];
	m1112 = Area_over_six*R3R3[2][3];
	m1210 = Area_over_six*R3R3[3][1];
	m1211 = Area_over_six*R3R3[3][2];
	
	//Bloco M34
	double m1013, m1014, m1015, m1016, m1113, m1114, m1115, m1116;

	m1013 = Area_3_over_twenty*R3[1][0];	
	m1014 = Area_3_over_twenty*R3[1][1];	
	m1015 = Area_3_over_twenty*R3[1][2];	
	m1016 = Area_3_over_twenty*R3[1][3];	
	m1113 = Area_3_over_twenty*R3[2][0];	
	m1114 = Area_3_over_twenty*R3[2][1];	
	m1115 = Area_3_over_twenty*R3[2][2];	
	m1116 = Area_3_over_twenty*R3[2][3];	

	//Reajustando a matriz M1 = Mhh + alpha*dt*Khh por blocos:

	//	|								|
	//	| M11 + alpha*dt*K11   M12 + alpha*dt*K12   M13 + alpha*dt*K13  |
	//   	| 								|
	//  M1 =| M12 + alpha*dt*K21   M22 + alpha*dt*K22   M23 + alpha*dt*K23	|
	//	|								|	
	//	| M13 + alpha*dt*K31   M23 + alpha*dt*K32   M33 + alpha*dt*K33	|
	//	|								|


	
	// Bloco 11
	M1[0][1] = m12 + alpha_dt*k12;
	M1[0][2] = m13 + alpha_dt*k13;
	M1[1][0] = m21 + alpha_dt*k21;
	M1[1][1] = m22 + alpha_dt*k22;
	M1[1][2] = m23 + alpha_dt*k23;
	M1[1][3] = m24 + alpha_dt*k24;
	M1[2][0] = m31 + alpha_dt*k31;
	M1[2][1] = m32 + alpha_dt*k32;
	M1[2][2] = m33 + alpha_dt*k33;
	M1[2][3] = m34 + alpha_dt*k34;
	M1[3][1] = m42 + alpha_dt*k42;
	M1[3][2] = m43 + alpha_dt*k43;

	//Bloco 12
	M1[0][5] = m16 + alpha_dt*k16;
	M1[0][6] = m17 + alpha_dt*k17;
	M1[1][4] = m25 + alpha_dt*k25;
	M1[1][5] = m26 + alpha_dt*k26;
	M1[1][6] = m27 + alpha_dt*k27;
	M1[1][7] = m28 + alpha_dt*k28;
	M1[2][4] = m35 + alpha_dt*k35;
	M1[2][5] = m36 + alpha_dt*k36;
	M1[2][6] = m37 + alpha_dt*k37;
	M1[2][7] = m38 + alpha_dt*k38;
	M1[3][5] = m46 + alpha_dt*k46;
	M1[3][6] = m47 + alpha_dt*k47;

	//Bloco 13
	M1[0][9 ] = m110 + alpha_dt*k110;
	M1[0][10] = m111 + alpha_dt*k1_11;
	M1[1][8 ] = m29  + alpha_dt*k29;
	M1[1][9 ] = m210 + alpha_dt*k210;
	M1[1][10] = m211 + alpha_dt*k211;
	M1[1][11] = m212 + alpha_dt*k212;
	M1[2][8 ] = m39  + alpha_dt*k39;
	M1[2][9 ] = m310 + alpha_dt*k310;
	M1[2][10] = m311 + alpha_dt*k311;
	M1[2][11] = m312 + alpha_dt*k312;
	M1[3][9 ] = m410 + alpha_dt*k410;
	M1[3][10] = m411 + alpha_dt*k411;

	//Bloco 21	
	M1[4][1] = m16 + alpha_dt*k52;
	M1[4][2] = m17 + alpha_dt*k53;
	M1[5][0] = m25 + alpha_dt*k61;
	M1[5][1] = m26 + alpha_dt*k62;
	M1[5][2] = m27 + alpha_dt*k63;
	M1[5][3] = m28 + alpha_dt*k64;
	M1[6][0] = m35 + alpha_dt*k71;
	M1[6][1] = m36 + alpha_dt*k72;
	M1[6][2] = m37 + alpha_dt*k73;
	M1[6][3] = m38 + alpha_dt*k74;
	M1[7][1] = m46 + alpha_dt*k82;
	M1[7][2] = m47 + alpha_dt*k83;

	//Bloco 22
	M1[4][5] = m56 + alpha_dt*k56;
	M1[4][6] = m57 + alpha_dt*k57;
	M1[5][4] = m65 + alpha_dt*k65;
	M1[5][5] = m66 + alpha_dt*k66;
	M1[5][6] = m67 + alpha_dt*k67;
	M1[5][7] = m68 + alpha_dt*k68;
	M1[6][4] = m75 + alpha_dt*k75;
	M1[6][5] = m76 + alpha_dt*k76;
	M1[6][6] = m77 + alpha_dt*k77;
	M1[6][7] = m78 + alpha_dt*k78;
	M1[7][5] = m86 + alpha_dt*k86;
	M1[7][6] = m87 + alpha_dt*k87;

	//Bloco 23
	M1[4][9 ] = m510 + alpha_dt*k510;
	M1[4][10] = m511 + alpha_dt*k511;
	M1[5][8 ] = m69  + alpha_dt*k69;
	M1[5][9 ] = m610 + alpha_dt*k610;
	M1[5][10] = m611 + alpha_dt*k611;
	M1[5][11] = m612 + alpha_dt*k612;
	M1[6][8 ] = m79  + alpha_dt*k79;
	M1[6][9 ] = m710 + alpha_dt*k710;
	M1[6][10] = m711 + alpha_dt*k711;
	M1[6][11] = m712 + alpha_dt*k712;
	M1[7][9 ] = m810 + alpha_dt*k810;
	M1[7][10] = m811 + alpha_dt*k811;

	//Bloco 31
	M1[8 ][1] = m110 + alpha_dt*k92;
	M1[8 ][2] = m111 + alpha_dt*k93;
	M1[9 ][0] = m29  + alpha_dt*k101;
	M1[9 ][1] = m210 + alpha_dt*k102;
	M1[9 ][2] = m211 + alpha_dt*k103;
	M1[9 ][3] = m212 + alpha_dt*k104;
	M1[10][0] = m39  + alpha_dt*k11_1;
	M1[10][1] = m310 + alpha_dt*k11_2;
	M1[10][2] = m311 + alpha_dt*k113;
	M1[10][3] = m312 + alpha_dt*k114;
	M1[11][1] = m410 + alpha_dt*k122;
	M1[11][2] = m411 + alpha_dt*k123;

	//Bloco 32
	M1[8 ][5] = m510 + alpha_dt*k96;
	M1[8 ][6] = m511 + alpha_dt*k97;
	M1[9 ][4] = m69  + alpha_dt*k105;
	M1[9 ][5] = m610 + alpha_dt*k106;
	M1[9 ][6] = m611 + alpha_dt*k107;
	M1[9 ][7] = m612 + alpha_dt*k108;
	M1[10][4] = m79  + alpha_dt*k115;
	M1[10][5] = m710 + alpha_dt*k116;
	M1[10][6] = m711 + alpha_dt*k117;
	M1[10][7] = m712 + alpha_dt*k118;
	M1[11][5] = m810 + alpha_dt*k126;
	M1[11][6] = m811 + alpha_dt*k127;

	//Bloco 33
	M1[8 ][9 ] = m910  + alpha_dt*k910; 
	M1[8 ][10] = m911  + alpha_dt*k911; 
	M1[9 ][8 ] = m109  + alpha_dt*k109;
	M1[9 ][9 ] = m1010 + alpha_dt*k1010;
	M1[9 ][10] = m1011 + alpha_dt*k1011;
	M1[9 ][11] = m1012 + alpha_dt*k1012;
	M1[10][8 ] = m119  + alpha_dt*k119;
	M1[10][9 ] = m1110 + alpha_dt*k1110;
	M1[10][10] = m1111 + alpha_dt*k1111;
	M1[10][11] = m1112 + alpha_dt*k1112;
	M1[11][9 ] = m1210 + alpha_dt*k1210;
	M1[11][10] = m1211 + alpha_dt*k1211;

	//Reajustando a matriz N1 = MhB + alpha*dt*KhB por blocos
			
	//	|		     |
	//	| M14 + alpha*dt*K14 |
	//  	| 		     |
	// N1 =	| M24 + alpha*dt*K24 |
	//	|		     |	
	//	| M34 + alpha*dt*K34 |	
	//	|		     |	

	//Bloco 14
	N1[1][0] = m213 + alpha_dt*k213;
	N1[1][1] = m214 + alpha_dt*k214;
	N1[1][2] = m215 + alpha_dt*k215;
	N1[1][3] = m216 + alpha_dt*k216;
	N1[2][0] = m313 + alpha_dt*k313;
	N1[2][1] = m314 + alpha_dt*k314;
	N1[2][2] = m315 + alpha_dt*k315;
	N1[2][3] = m316 + alpha_dt*k316;
	
	//Bloco 24
	N1[5][0] = m613 + alpha_dt*k613;
	N1[5][1] = m614 + alpha_dt*k614;
	N1[5][2] = m615 + alpha_dt*k615;
	N1[5][3] = m616 + alpha_dt*k616;
	N1[6][0] = m713 + alpha_dt*k713;
	N1[6][1] = m714 + alpha_dt*k714;
	N1[6][2] = m715 + alpha_dt*k715;
	N1[6][3] = m716 + alpha_dt*k716;

	//Bloco 34
	N1[9 ][0] = m1013 + alpha_dt*k1013;
	N1[9 ][1] = m1014 + alpha_dt*k1014;
	N1[9 ][2] = m1015 + alpha_dt*k1015;
	N1[9 ][3] = m1016 + alpha_dt*k1016;
	N1[10][0] = m1113 + alpha_dt*k1113;
	N1[10][1] = m1114 + alpha_dt*k1114;
	N1[10][2] = m1115 + alpha_dt*k1115;
	N1[10][3] = m1116 + alpha_dt*k1116;

	//Reajustando a matriz M2 = MBh + alpha*dt*KBh por blocos
	//  	| 		    			   			 |
	// M2 =	| M14 + alpha*dt*K41	M24 + alpha*dt*K42    M34 + alpha*dt*K43 |
	//	|		    						 |	
	
	//Bloco 41
	M2[0][1] = 	      alpha_dt*k132;
	M2[1][1] = m214 + alpha_dt*k142;
	M2[2][1] = m314 + alpha_dt*k152;
	M2[3][1] = 	      alpha_dt*k162;
	M2[0][2] = 	      alpha_dt*k133;
	M2[1][2] = m215 + alpha_dt*k143;
	M2[2][2] = m315 + alpha_dt*k153;
	M2[3][2] = 	  alpha_dt*k163;

	//Bloco 42
	M2[0][5] = 	      alpha_dt*k136;
	M2[1][5] = m614 + alpha_dt*k146;
	M2[2][5] = m714 + alpha_dt*k156;
	M2[3][5] = 	      alpha_dt*k166;
	M2[0][6] =        alpha_dt*k137;
	M2[1][6] = m615 + alpha_dt*k147;
	M2[2][6] = m715 + alpha_dt*k157;
	M2[3][6] = 	      alpha_dt*k167;
	
	//Bloco 43
	M2[0][9]  = 	     alpha_dt*k1310;
	M2[1][9]  = m1014 + alpha_dt*k1410;
	M2[2][9]  = m1114 + alpha_dt*k1510;
	M2[3][9]  = 	     alpha_dt*k1610;
	M2[0][10] =         alpha_dt*k1311;
	M2[1][10] = m1015 + alpha_dt*k1411;
	M2[2][10] = m1115 + alpha_dt*k1511;
	M2[3][10] = 	     alpha_dt*k1611;
	
	//Calculo do Res1
	double MhB = Area_3_over_twenty;

	double KhB11 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[0][0] + (x32 + 2*x13 + 2*x21)*R1Ay[0][0]);		
	double KhB12 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[0][1] + (x32 + 2*x13 + 2*x21)*R1Ay[0][1]);		
	double KhB13 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[0][2] + (x32 + 2*x13 + 2*x21)*R1Ay[0][2]);		
	double KhB14 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[0][3] + (x32 + 2*x13 + 2*x21)*R1Ay[0][3]);		
	double KhB41 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[3][0] + (x32 + 2*x13 + 2*x21)*R1Ay[3][0]);		
	double KhB42 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[3][1] + (x32 + 2*x13 + 2*x21)*R1Ay[3][1]);		
	double KhB43 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[3][2] + (x32 + 2*x13 + 2*x21)*R1Ay[3][2]);		
	double KhB44 = nine_over_forty*((y23 + 2*y31 + 2*y12)*R1Ax[3][3] + (x32 + 2*x13 + 2*x21)*R1Ay[3][3]);		
	
	double KhB51 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[0][0] + (2*x32 + x13 + 2*x21)*R2Ay[0][0]);		
	double KhB52 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[0][1] + (2*x32 + x13 + 2*x21)*R2Ay[0][1]);		
	double KhB53 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[0][2] + (2*x32 + x13 + 2*x21)*R2Ay[0][2]);		
	double KhB54 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[0][3] + (2*x32 + x13 + 2*x21)*R2Ay[0][3]);		
	double KhB81 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[3][0] + (2*x32 + x13 + 2*x21)*R2Ay[3][0]);		
	double KhB82 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[3][1] + (2*x32 + x13 + 2*x21)*R2Ay[3][1]);		
	double KhB83 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[3][2] + (2*x32 + x13 + 2*x21)*R2Ay[3][2]);		
	double KhB84 = nine_over_forty*((2*y23 + y31 + 2*y12)*R2Ax[3][3] + (2*x32 + x13 + 2*x21)*R2Ay[3][3]);		

	double KhB91  = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[0][0] + (2*x32 + 2*x13 + x21)*R3Ay[0][0]);		
	double KhB92  = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[0][1] + (2*x32 + 2*x13 + x21)*R3Ay[0][1]);		
	double KhB93  = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[0][2] + (2*x32 + 2*x13 + x21)*R3Ay[0][2]);		
	double KhB94  = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[0][3] + (2*x32 + 2*x13 + x21)*R3Ay[0][3]);		
	double KhB121 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[3][0] + (2*x32 + 2*x13 + x21)*R3Ay[3][0]);		
	double KhB122 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[3][0] + (2*x32 + 2*x13 + x21)*R3Ay[3][0]);		
	double KhB123 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[3][0] + (2*x32 + 2*x13 + x21)*R3Ay[3][0]);		
	double KhB124 = nine_over_forty*((2*y23 + 2*y31 + y12)*R3Ax[3][0] + (2*x32 + 2*x13 + x21)*R3Ay[3][0]);		

	Res1[0]  = - MhB*dUeB[0]                   - (KhB11*UeB[0]  + KhB12*UeB[1]  + KhB13*UeB[2]  + KhB14*UeB[3]);
	Res1[1]  = - (m214*dUeB[1] + m215*dUeB[2]) - (k213*UeB[0]   + k214*UeB[1]   + k215*UeB[2]   + k216*UeB[3] );
	Res1[2]  = - (m314*dUeB[1] + m315*dUeB[2]) - (k313*UeB[0]   + k314*UeB[1]   + k315*UeB[2]   + k316*UeB[3] );
	Res1[3]  = - MhB*dUeB[3]                   - (KhB41*UeB[0]  + KhB42*UeB[1]  + KhB43*UeB[2]  + KhB44*UeB[3]);

	Res1[4]  = - MhB*dUeB[0]                   - (KhB51*UeB[0]  + KhB52*UeB[1]  + KhB53*UeB[2]  + KhB54*UeB[3]);
	Res1[5]  = - (m614*dUeB[1] + m615*dUeB[2]) - (k613*UeB[0]   + k614*UeB[1]   + k615*UeB[2]   + k615*UeB[3] );
	Res1[6]  = - (m714*dUeB[1] + m715*dUeB[2]) - (k713*UeB[0]   + k714*UeB[1]   + k715*UeB[2]   + k716*UeB[3] );
	Res1[7]  = - MhB*dUeB[3]                   - (KhB81*UeB[0]  + KhB82*UeB[1]  + KhB83*UeB[2]  + KhB84*UeB[3]);

	Res1[8]  = - MhB*dUeB[0]                     - (KhB91*UeB[0]  + KhB92*UeB[1]  + KhB93*UeB[2]  + KhB94*UeB[3] );
	Res1[9]  = - (m1014*dUeB[1] + m1015*dUeB[2]) - (k1013*UeB[0]  + k1014*UeB[1]  + k1015*UeB[2]  + k1016*UeB[3] );
	Res1[10] = - (m1114*dUeB[1] + m1115*dUeB[2]) - (k1113*UeB[0]  + k1114*UeB[1]  + k1115*UeB[2]  + k1116*UeB[3] );
	Res1[11] = - MhB*dUeB[3]                     - (KhB121*UeB[0] + KhB122*UeB[1] + KhB123*UeB[2] + KhB124*UeB[3]);
	
 	double k11   = one_sixth*(y23*R1AxR1[0][0] + x32*R1AyR1[0][0]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[0][0];
 	double k14   = one_sixth*(y23*R1AxR1[0][3] + x32*R1AyR1[0][3]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[0][3];
 	double k41   = one_sixth*(y23*R1AxR1[3][0] + x32*R1AyR1[3][0]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[3][0];
 	double k44   = one_sixth*(y23*R1AxR1[3][3] + x32*R1AyR1[3][3]) + delta_quarter_Area*(y23*y23 + x32*x32)*R1R1[3][3];

	double k15   = one_sixth*(y31*R1AxR2[0][0] + x13*R1AyR2[0][0]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[0][0];
	double k18   = one_sixth*(y31*R1AxR2[0][3] + x13*R1AyR2[0][3]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[0][3];
	double k45   = one_sixth*(y31*R1AxR2[3][0] + x13*R1AyR2[3][0]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[3][0];
	double k48   = one_sixth*(y31*R1AxR2[3][3] + x13*R1AyR2[3][3]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[3][3];

	double k19    = one_sixth*(y12*R1AxR3[0][0] + x21*R1AyR3[0][0]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[0][0];
	double k1_12  = one_sixth*(y12*R1AxR3[0][3] + x21*R1AyR3[0][3]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[0][3];
	double k49    = one_sixth*(y12*R1AxR3[3][0] + x21*R1AyR3[3][0]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[3][0];
	double k412   = one_sixth*(y12*R1AxR3[3][3] + x21*R1AyR3[3][3]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[3][3];

	double k51   = one_sixth*(y23*R2AxR1[0][0] + x32*R2AyR1[0][0]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[0][0];
	double k54   = one_sixth*(y23*R2AxR1[0][3] + x32*R2AyR1[0][3]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[0][3];
	double k81   = one_sixth*(y23*R2AxR1[3][0] + x32*R2AyR1[3][0]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[3][0];
	double k84   = one_sixth*(y23*R2AxR1[3][3] + x32*R2AyR1[3][3]) + delta_quarter_Area*(y23*y31 + x32*x13)*R1R2[3][3];

	double k55   = one_sixth*(y31*R2AxR2[0][0] + x13*R2AyR2[0][0]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[0][0];
	double k58   = one_sixth*(y31*R2AxR2[0][3] + x13*R2AyR2[0][3]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[0][3];
	double k85   = one_sixth*(y31*R2AxR2[3][0] + x13*R2AyR2[3][0]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[3][0];
	double k88   = one_sixth*(y31*R2AxR2[3][3] + x13*R2AyR2[3][3]) + delta_quarter_Area*(y31*y31 + x13*x13)*R2R2[3][3];

	double k59   = one_sixth*(y12*R2AxR3[0][0] + x21*R2AyR3[0][0]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[0][0];
	double k512  = one_sixth*(y12*R2AxR3[0][3] + x21*R2AyR3[0][3]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[0][3];
	double k89   = one_sixth*(y12*R2AxR3[3][0] + x21*R2AyR3[3][0]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[3][0];
	double k812  = one_sixth*(y12*R2AxR3[3][3] + x21*R2AyR3[3][3]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[3][3];

	double k91   = one_sixth*(y23*R3AxR1[0][0] + x32*R3AyR1[0][0]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[0][0];
	double k94   = one_sixth*(y23*R3AxR1[0][3] + x32*R3AyR1[0][3]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[0][3];
	double k121  = one_sixth*(y23*R3AxR1[3][0] + x32*R3AyR1[3][0]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[3][0];
	double k124  = one_sixth*(y23*R3AxR1[3][3] + x32*R3AyR1[3][3]) + delta_quarter_Area*(y23*y12 + x32*x21)*R1R3[3][3];

	double k95   = one_sixth*(y31*R3AxR2[0][0] + x13*R3AyR2[0][0]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[0][0];
	double k98   = one_sixth*(y31*R3AxR2[0][3] + x13*R3AyR2[0][3]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[0][3];
	double k125  = one_sixth*(y31*R3AxR2[3][0] + x13*R3AyR2[3][0]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[3][0];
	double k128  = one_sixth*(y31*R3AxR2[3][3] + x13*R3AyR2[3][3]) + delta_quarter_Area*(y31*y12 + x13*x21)*R2R3[3][3];

	double k99   = one_sixth*(y12*R3AxR3[0][0] + x21*R3AyR3[0][0]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[0][0];
	double k912  = one_sixth*(y12*R3AxR3[0][3] + x21*R3AyR3[0][3]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[0][3];
	double k129  = one_sixth*(y12*R3AxR3[3][0] + x21*R3AyR3[3][0]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[3][0];
	double k1212 = one_sixth*(y12*R3AxR3[3][3] + x21*R3AyR3[3][3]) + delta_quarter_Area*(y12*y12 + x21*x21)*R3R3[3][3];
	
	
	double KhhUe[12];
	
	KhhUe[0]  = k11*Ue[0] + k12*Ue[1] + k13*Ue[2] + k14*Ue[3] + k15*Ue[4] + k16*Ue[5] + k17*Ue[6] + k18*Ue[7] + k19*Ue[8] + k110*Ue[9] + k1_11*Ue[10] + k1_12*Ue[11];
	KhhUe[1]  = k21*Ue[0] + k22*Ue[1] + k23*Ue[2] + k24*Ue[3] + k25*Ue[4] + k26*Ue[5] + k27*Ue[6] + k28*Ue[7] + k29*Ue[8] + k210*Ue[9] + k211 *Ue[10] + k212 *Ue[11];
	KhhUe[2]  = k31*Ue[0] + k32*Ue[1] + k33*Ue[2] + k34*Ue[3] + k35*Ue[4] + k36*Ue[5] + k37*Ue[6] + k38*Ue[7] + k39*Ue[8] + k310*Ue[9] + k311 *Ue[10] + k312 *Ue[11];
	KhhUe[3]  = k41*Ue[0] + k42*Ue[1] + k43*Ue[2] + k44*Ue[3] + k45*Ue[4] + k46*Ue[5] + k47*Ue[6] + k48*Ue[7] + k49*Ue[8] + k410*Ue[9] + k411 *Ue[10] + k412 *Ue[11];

	KhhUe[4]  = k51*Ue[0] + k52*Ue[1] + k53*Ue[2] + k54*Ue[3] + k55*Ue[4] + k56*Ue[5] + k57*Ue[6] + k58*Ue[7] + k59*Ue[8] + k510*Ue[9] + k511* Ue[10] +  k512*Ue[11];
	KhhUe[5]  = k61*Ue[0] + k62*Ue[1] + k63*Ue[2] + k64*Ue[3] + k65*Ue[4] + k66*Ue[5] + k67*Ue[6] + k68*Ue[7] + k69*Ue[8] + k610*Ue[9] + k611* Ue[10] +  k612*Ue[11];
	KhhUe[6]  = k71*Ue[0] + k72*Ue[1] + k73*Ue[2] + k74*Ue[3] + k75*Ue[4] + k76*Ue[5] + k77*Ue[6] + k78*Ue[7] + k79*Ue[8] + k710*Ue[9] + k711* Ue[10] +  k712*Ue[11];
	KhhUe[7]  = k81*Ue[0] + k82*Ue[1] + k83*Ue[2] + k84*Ue[3] + k85*Ue[4] + k86*Ue[5] + k87*Ue[6] + k88*Ue[7] + k89*Ue[8] + k810*Ue[9] + k811* Ue[10] +  k812*Ue[11];

	KhhUe[8 ]  = k91 * Ue[0] + k92*Ue[1]   + k93*Ue[2]  + k94*Ue[3]  + k95*Ue[4]  + k96*Ue[5]  + k97*Ue[6]  + k98*Ue[7]  + k99*Ue[8]  + k910*Ue[9]  + k911* Ue[10]  +  k912*Ue[11];
	KhhUe[9 ]  = k101* Ue[0] + k102*Ue[1]  + k103*Ue[2] + k104*Ue[3] + k105*Ue[4] + k106*Ue[5] + k107*Ue[6] + k108*Ue[7] + k109*Ue[8] + k1010*Ue[9] + k1011* Ue[10] +  k1012*Ue[11];
	KhhUe[10]  = k11_1*Ue[0] + k11_2*Ue[1] + k113*Ue[2] + k114*Ue[3] + k115*Ue[4] + k116*Ue[5] + k117*Ue[6] + k118*Ue[7] + k119*Ue[8] + k1110*Ue[9] + k1111* Ue[10] +  k1112*Ue[11];
	KhhUe[11]  = k121* Ue[0] + k122*Ue[1]  + k123*Ue[2] + k124*Ue[3] + k125*Ue[4] + k126*Ue[5] + k127*Ue[6] + k128*Ue[7] + k129*Ue[8] + k1210*Ue[9] + k1211* Ue[10] +  k1212*Ue[11];

	
	double Mhh1 = Area_over_six;
	double Mhh2 = Area_over_twelve;
	
	Res1[0]  += - (Mhh1*dUe[0] + Mhh2*dUe[4] + Mhh2*dUe[8])                                          - KhhUe[0];
	Res1[1]  += - (m22*dUe[1]  + m23*dUe[2] + m26*dUe[5] + m27*dUe[6] + m210*dUe[9 ] + m211*dUe[10]) - KhhUe[1];
	Res1[2]  += - (m32*dUe[1]  + m33*dUe[2] + m36*dUe[5] + m37*dUe[6] + m310*dUe[10] + m311*dUe[10]) - KhhUe[2];
	Res1[3]  += - (Mhh1*dUe[3] + Mhh2*dUe[7] + Mhh2*dUe[11]) 					                         - KhhUe[3];

	Res1[4]  += - (Mhh2*dUe[0] + Mhh1*dUe[4] + Mhh2*dUe[8])                                          - KhhUe[4];
	Res1[5]  += - (m26*dUe[1]  + m27*dUe[2] + m66*dUe[5] + m67*dUe[6] + m610*dUe[9 ] + m611*dUe[10]) - KhhUe[5];
	Res1[6]  += - (m36*dUe[1]  + m37*dUe[2] + m76*dUe[5] + m77*dUe[6] + m710*dUe[10] + m711*dUe[10]) - KhhUe[6];
	Res1[7]  += - (Mhh2*dUe[3] + Mhh1*dUe[7] + Mhh2*dUe[11]) 					                         - KhhUe[7];

	Res1[8]   += - (Mhh2*dUe[0]  + Mhh2*dUe[4] + Mhh1*dUe[8])                                               - KhhUe[8];
	Res1[9]   += - (m210*dUe[1]  + m211*dUe[2] + m610*dUe[5] + m611*dUe[6] + m1010*dUe[9 ] + m1011*dUe[10]) - KhhUe[9];
	Res1[10]  += - (m310*dUe[1]  + m311*dUe[2] + m710*dUe[5] + m711*dUe[6] + m1110*dUe[10] + m1111*dUe[10]) - KhhUe[10];
	Res1[11]  += - (Mhh2*dUe[3]  + Mhh2*dUe[7] + Mhh1*dUe[11])                            	                 - KhhUe[11];

	//Calculo Res2
	double KBB = (81.0*delta*(y23*y23 + y31*y31 + y12*y12 + y23*y31 + y23*y12 + y31*y12 + x32*x32 + x13*x13 + x21*x21 + x32*x13 + x32*x21 + x13*x21))/(40.0*Area);
	double MBB = (81.0 * Area) / 280.0;

	Res2[0]  = - MBB*dUeB[0] - KBB*UeB[0];
	Res2[1]  = - MBB*dUeB[1] - KBB*UeB[1];
	Res2[2]  = - MBB*dUeB[2] - KBB*UeB[2];
	Res2[3]  = - MBB*dUeB[3] - KBB*UeB[3];
	
	double KBhUe[4];
	
	double KBh11 = nine_over_forty*(y23*AxR1[0][0] + x32*AyR1[0][0]);
	double KBh21 = nine_over_forty*(y23*AxR1[1][0] + x32*AyR1[1][0]);
	double KBh31 = nine_over_forty*(y23*AxR1[2][0] + x32*AyR1[2][0]);
	double KBh41 = nine_over_forty*(y23*AxR1[3][0] + x32*AyR1[3][0]);
	double KBh14 = nine_over_forty*(y23*AxR1[0][3] + x32*AyR1[0][3]);
	double KBh24 = nine_over_forty*(y23*AxR1[1][3] + x32*AyR1[1][3]);
	double KBh34 = nine_over_forty*(y23*AxR1[2][3] + x32*AyR1[2][3]);
	double KBh44 = nine_over_forty*(y23*AxR1[3][3] + x32*AyR1[3][3]);

	double KBh15 = nine_over_forty*(y31*AxR2[0][0] + x13*AyR2[0][0]);
	double KBh25 = nine_over_forty*(y31*AxR2[1][0] + x13*AyR2[1][0]);
	double KBh35 = nine_over_forty*(y31*AxR2[2][0] + x13*AyR2[2][0]);
	double KBh45 = nine_over_forty*(y31*AxR2[3][0] + x13*AyR2[3][0]);
	double KBh18 = nine_over_forty*(y31*AxR2[0][3] + x13*AyR2[0][3]);
	double KBh28 = nine_over_forty*(y31*AxR2[1][3] + x13*AyR2[1][3]);
	double KBh38 = nine_over_forty*(y31*AxR2[2][3] + x13*AyR2[2][3]);
	double KBh48 = nine_over_forty*(y31*AxR2[3][3] + x13*AyR2[3][3]);

	double KBh19  = nine_over_forty*(y12*AxR3[0][0] + x21*AyR3[0][0]);
	double KBh29  = nine_over_forty*(y12*AxR3[1][0] + x21*AyR3[1][0]);
	double KBh39  = nine_over_forty*(y12*AxR3[2][0] + x21*AyR3[2][0]);
	double KBh49  = nine_over_forty*(y12*AxR3[3][0] + x21*AyR3[3][0]);
	double KBh112 = nine_over_forty*(y12*AxR3[0][3] + x21*AyR3[0][3]);
	double KBh212 = nine_over_forty*(y12*AxR3[1][3] + x21*AyR3[1][3]);
	double KBh312 = nine_over_forty*(y12*AxR3[2][3] + x21*AyR3[2][3]);
	double KBh412 = nine_over_forty*(y12*AxR3[3][3] + x21*AyR3[3][3]);
	
	KBhUe[0]  = KBh11*Ue[0] + k132*Ue[1] + k133*Ue[2] + KBh14*Ue[3] + KBh15*Ue[4] + k136*Ue[5] + k137*Ue[6] + KBh18*Ue[7] + KBh19*Ue[8] + k1310*Ue[9] + k1311*Ue[10] + KBh112*Ue[11];
	KBhUe[1]  = KBh21*Ue[0] + k142*Ue[1] + k143*Ue[2] + KBh24*Ue[3] + KBh25*Ue[4] + k146*Ue[5] + k147*Ue[6] + KBh28*Ue[7] + KBh29*Ue[8] + k1410*Ue[9] + k1411*Ue[10] + KBh212*Ue[11];
	KBhUe[2]  = KBh31*Ue[0] + k152*Ue[1] + k153*Ue[2] + KBh34*Ue[3] + KBh35*Ue[4] + k156*Ue[5] + k157*Ue[6] + KBh38*Ue[7] + KBh39*Ue[8] + k1510*Ue[9] + k1511*Ue[10] + KBh312*Ue[11];
	KBhUe[3]  = KBh41*Ue[0] + k162*Ue[1] + k163*Ue[2] + KBh44*Ue[3] + KBh45*Ue[4] + k166*Ue[5] + k167*Ue[6] + KBh48*Ue[7] + KBh49*Ue[8] + k1610*Ue[9] + k1611*Ue[10] + KBh412*Ue[11];

	Res2[0]  += - (MhB*dUe[0] + MhB *dUe[4]  + MhB*dUe[8])- KBhUe[0];
	Res2[1]  += - (m214*dUe[1] + m215*dUe[2] + m614*dUe[5] + m615*dUe[6] + m1014*dUe[9] + m1015*dUe[10])  - KBhUe[1];
	Res2[2]  += - (m314*dUe[1] + m315*dUe[2] + m714*dUe[5] + m715*dUe[6] + m1114*dUe[9] + m1115*dUe[10])  - KBhUe[2];
	Res2[3]  += - (MhB*dUe[3] + MhB*dUe[7] + MhB*dUe[11]) - KBhUe[3];
	
}
