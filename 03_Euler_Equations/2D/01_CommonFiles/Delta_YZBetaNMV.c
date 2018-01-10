#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

double Delta_YZBetaNMV(double tolerance, double *delta_old_NMV, double *gradUx, double *gradUy, double (*Ax)[4], double (*Ay)[4], double (*A0)[4], 
				double *dUb, double y23, double y31, double y12, double x32, double x13, double x21, double twoArea, int e, double *invY, double *Ub)
{
	double beta, *AxgradUx, *AygradUy, *Z, norm_YZ, norm_YdUx, norm_YdUy, sum_norm_YdU, norm_YU, delta1, delta;
	double J1, J2, norm_gradRHO, Hyzbeta, aux_sum_norm_YdU, aux_Hyzbeta, aux_norm_YU;
	int i;

	AxgradUx = (double*) mycalloc("AxgradUx of 'Delta_YZBeta'", 4, sizeof(double));
	AygradUy = (double*) mycalloc("AygradUy of 'Delta_YZBeta'", 4, sizeof(double));
	Z = (double*) mycalloc("Z of 'Delta_YZBeta'", 4, sizeof(double));
	
	// *** Ax*gradUx and Ay*gradUy coefficients
	AxgradUx[0] = gradUx[1];
	AxgradUx[1] = Ax[1][0]*gradUx[0] + Ax[1][1]*gradUx[1] + Ax[1][2]*gradUx[2] + Ax[1][3]*gradUx[3];
	AxgradUx[2] = Ax[2][0]*gradUx[0] + Ax[2][1]*gradUx[1] + Ax[2][2]*gradUx[2];
	AxgradUx[3] = Ax[3][0]*gradUx[0] + Ax[3][1]*gradUx[1] + Ax[3][2]*gradUx[2] + Ax[3][3]*gradUx[3];

	AygradUy[0] = gradUy[2];
	AygradUy[1] = Ay[1][0]*gradUy[0] + Ay[1][1]*gradUy[1] + Ay[1][2]*gradUy[2];
	AygradUy[2] = Ay[2][0]*gradUy[0] + Ay[2][1]*gradUy[1] + Ay[2][2]*gradUy[2] + Ay[2][3]*gradUy[3];
	AygradUy[3] = Ay[3][0]*gradUy[0] + Ay[3][1]*gradUy[1] + Ay[3][2]*gradUy[2] + Ay[3][3]*gradUy[3];
/*	#ifdef debug
		for(i = 0; i < 4; i++){
			printf("AxgradUx[%d] = %lf \t AygradUy[%d] = %lf \n", i, AxgradUx[i], i, AygradUy[i]);
		}
		getchar();
	#endif */
	
	// Calculo de Z = dU + AxgradUx + AygradUy
	for (i = 0; i < 4; i++)
		Z[i] = dUb[i] + AxgradUx[i] + AygradUy[i];
	
	// norm_YZ = ||Y^{-1}*Z||
	norm_YZ = sqrt(invY[0]*Z[0]*invY[0]*Z[0] + invY[1]*Z[1]*invY[1]*Z[1] + invY[2]*Z[2]*invY[2]*Z[2] + invY[3]*Z[3]*invY[3]*Z[3]);
	
	// norm_YdUx = ||Y^{-1}*dU/dx||^2
	norm_YdUx = invY[0]*gradUx[0]*invY[0]*gradUx[0] + invY[1]*gradUx[1]*invY[1]*gradUx[1] + invY[2]*gradUx[2]*invY[2]*gradUx[2] + invY[3]*gradUx[3]*invY[3]*gradUx[3];
	
	//norm_YdUy = ||Y^{-1}*dU/dy||^2
	norm_YdUy = invY[0]*gradUy[0]*invY[0]*gradUy[0] + invY[1]*gradUy[1]*invY[1]*gradUy[1] + invY[2]*gradUy[2]*invY[2]*gradUy[2] + invY[3]*gradUy[3]*invY[3]*gradUy[3];
	
	// sum_norm_YdU = ||Y^{-1}*dU/dx||^2 + ||Y^{-1}*dU/dy||^2
	sum_norm_YdU = norm_YdUx + norm_YdUy;
	
	// norm_YU = ||Y^{-1}*U||
	norm_YU = sqrt(invY[0]*Ub[0]*invY[0]*Ub[0] + invY[1]*Ub[1]*invY[1]*Ub[1] + invY[2]*Ub[2]*invY[2]*Ub[2] + invY[3]*Ub[3]*invY[3]*Ub[3]);
	
	// J = gradRHO/||gradRHO||
	norm_gradRHO = sqrt(gradUx[0]*gradUx[0] + gradUy[0]*gradUy[0]);
	if(fabs(norm_gradRHO) >= 1e-12){
		J1 = gradUx[0]/norm_gradRHO;
		J2 = gradUy[0]/norm_gradRHO;
	}else{
		J1 = J2 = 0.0;
	}
	
	Hyzbeta = (fabs(J1*y23 + J2*x32) + fabs(J1*y31 + J2*x13) + fabs(J1*y12 + J2*x21))/twoArea;
	if(fabs(Hyzbeta) >= 1e-12){
		Hyzbeta = 2.0*(1.0/Hyzbeta);
	}else{
		Hyzbeta = 0.0;
	}

	//Calculando delta1
	beta = 1.0;
	if(fabs(sum_norm_YdU) >= 1e-12)
		aux_sum_norm_YdU = 1.0 / sqrt(sum_norm_YdU); 
	else		
		aux_sum_norm_YdU = 0.0;

	if(fabs(norm_YU) >= 1e-12)
		aux_norm_YU = 1.0; 
	else
		aux_norm_YU = 0.0; 

	aux_Hyzbeta = pow((Hyzbeta*0.5),beta);
	delta1 = norm_YZ*aux_sum_norm_YdU*aux_norm_YU*aux_Hyzbeta;

	double aux_delta, w = 0.5;	
	aux_delta = 0.5*delta1; 
	delta = w*(aux_delta) + (1 - w)*delta_old_NMV[e];
	delta_old_NMV[e] = delta;
	free(AxgradUx);
	free(AygradUy);
	free(Z);
	
	return delta;
		
}
