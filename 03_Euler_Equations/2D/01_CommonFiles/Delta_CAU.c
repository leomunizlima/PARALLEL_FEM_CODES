#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

double Delta_CAU(double tolerance, double *delta_old, double *gradUx, double *gradUy, double (*Ax)[4], double (*Ay)[4], double (*A0)[4],
				double *dUb, double y23, double y31, double y12, double x32, double x13, double x21, double twoArea, int e, double *invY, double *Ub)
{
	
	double *AxgradUx, *AygradUy, *gradksiUx, *gradksiUy, *LUh, delta, norma_gradU, norma_LUh, norma_gradksiU;
	int i;
	
	AxgradUx = (double*) mycalloc("AxgradUx of 'Delta_CAU'", 4, sizeof(double));
	AygradUy = (double*) mycalloc("AygradUy of 'Delta_CAU'", 4, sizeof(double));
	gradksiUx = (double*) mycalloc("gradksiUx of 'Delta_CAU'", 4, sizeof(double));
	gradksiUy = (double*) mycalloc("gradksiUy of 'Delta_CAU'", 4, sizeof(double));
	LUh = (double*) mycalloc("LUh of 'Delta_CAU'", 4, sizeof(double));
	
	// *** Coeficientes de Ax gradUx e Ay gradUy
	AxgradUx[0] = gradUx[1];
	AxgradUx[1] = Ax[1][0]*gradUx[0] + Ax[1][1]*gradUx[1] + Ax[1][2]*gradUx[2] + Ax[1][3]*gradUx[3];
	AxgradUx[2] = Ax[2][0]*gradUx[0] + Ax[2][1]*gradUx[1] + Ax[2][2]*gradUx[2];
	AxgradUx[3] = Ax[3][0]*gradUx[0] + Ax[3][1]*gradUx[1] + Ax[3][2]*gradUx[2] + Ax[3][3]*gradUx[3];

	AygradUy[0] = gradUy[2];
	AygradUy[1] = Ay[1][0]*gradUy[0] + Ay[1][1]*gradUy[1] + Ay[1][2]*gradUy[2];
	AygradUy[2] = Ay[2][0]*gradUy[0] + Ay[2][1]*gradUy[1] + Ay[2][2]*gradUy[2] + Ay[2][3]*gradUy[3];
	AygradUy[3] = Ay[3][0]*gradUy[0] + Ay[3][1]*gradUy[1] + Ay[3][2]*gradUy[2] + Ay[3][3]*gradUy[3];

	// *** Coeficientes de gradksiUx = J1i gradUxi e gradksiUy = J2i gradxi  
	for (i = 0; i < 4; i++)
	{
		gradksiUx[i] = (y23*gradUx[i] + x32*gradUy[i]) / twoArea;
		gradksiUy[i] = (y31*gradUx[i] + x13*gradUy[i]) / twoArea;   
	}

	// Norma de |grad ksi * Uh| em A0(-1)
	norma_gradksiU = (gradksiUx[0] * (A0[0][0]*gradksiUx[0] + 2.0*A0[0][1]*gradksiUx[1] + 2.0*A0[0][2]*gradksiUx[2] + 2.0*A0[0][3]*gradksiUx[3]) +
							gradksiUx[1] * (A0[1][1]*gradksiUx[1] + 2.0*A0[1][2]*gradksiUx[2] + 2.0*A0[1][3]*gradksiUx[3]) +
							gradksiUx[2] * (A0[2][2]*gradksiUx[2] + 2.0*A0[2][3]*gradksiUx[3]) + 
							gradksiUx[3] * A0[3][3] * gradksiUx[3]) +
							(gradksiUy[0] * (A0[0][0]*gradksiUy[0] + 2.0*A0[0][1]*gradksiUy[1] + 2.0*A0[0][2]*gradksiUy[2] + 2.0*A0[0][3]*gradksiUy[3]) +
							gradksiUy[1] * (A0[1][1]*gradksiUy[1] + 2.0*A0[1][2]*gradksiUy[2] + 2.0*A0[1][3]*gradksiUy[3]) +
							gradksiUy[2] * (A0[2][2]*gradksiUy[2] + 2.0*A0[2][3]*gradksiUy[3]) + 
							gradksiUy[3] * A0[3][3] * gradksiUy[3]); 

	norma_gradksiU = sqrt(norma_gradksiU);

	// Calculo de |LUh| = | dU + AxgradUx + AygradUy | em A0(-1)
	for (i = 0; i < 4; i++)
		LUh[i] = dUb[i] + AxgradUx[i] + AygradUy[i];

	norma_LUh = (LUh[0] * (A0[0][0]*LUh[0] + 2.0*A0[0][1]*LUh[1] + 2.0*A0[0][2]*LUh[2] + 2.0*A0[0][3]*LUh[3]) +
						LUh[1] * (A0[1][1]*LUh[1] + 2.0*A0[1][2]*LUh[2] + 2.0*A0[1][3]*LUh[3]) +
						LUh[2] * (A0[2][2]*LUh[2] + 2.0*A0[2][3]*LUh[3]) + 
						LUh[3] * A0[3][3] * LUh[3]);

	norma_LUh = sqrt(norma_LUh);

	// Norma de |grad Uh| em A0(-1)
	norma_gradU = (gradUx[0] * (A0[0][0]*gradUx[0] + 2.0*A0[0][1]*gradUx[1] + 2.0*A0[0][2]*gradUx[2] + 2.0*A0[0][3]*gradUx[3]) +
						gradUx[1] * (A0[1][1]*gradUx[1] + 2.0*A0[1][2]*gradUx[2] + 2.0*A0[1][3]*gradUx[3]) +
						gradUx[2] * (A0[2][2]*gradUx[2] + 2.0*A0[2][3]*gradUx[3]) + 
						gradUx[3] * A0[3][3] * gradUx[3]) +
						(gradUy[0] * (A0[0][0]*gradUy[0] + 2.0*A0[0][1]*gradUy[1] + 2.0*A0[0][2]*gradUy[2] + 2.0*A0[0][3]*gradUy[3]) +
						gradUy[1] * (A0[1][1]*gradUy[1] + 2.0*A0[1][2]*gradUy[2] + 2.0*A0[1][3]*gradUy[3]) +
						gradUy[2] * (A0[2][2]*gradUy[2] + 2.0*A0[2][3]*gradUy[3]) + 
						gradUy[3] * A0[3][3]*gradUy[3]);
						
	norma_gradU = sqrt(norma_gradU);

	// ... Operador de Captura  de Choque  CAU
	if (norma_gradU >= tolerance)
		delta = norma_LUh / norma_gradksiU;
	else
		delta = 0.0;
	
	free(AxgradUx);
	free(AygradUy);
	free(gradksiUx);
	free(gradksiUy);
	free(LUh);
	
	return delta;
		
}
