#include "TranspEquation.h"

double YZBeta_ShockCapture(double kx, double ky, double Be_x, double Be_y, double gamma, double ue1, double ue2, double ue3, double ueb, double dueb, double feb, 
                        double y23, double y31, double y12, double x32, double x13, double x21, double invArea, double h)
{
	double Eu, BetaGradu, Ru, normRu, Gradu[2], normGradu, normu;	

	BetaGradu = 0.5*invArea*(ue1*(Be_x*y23 + Be_y*x32) + ue2*(Be_x*y31 + Be_y*x13) + ue3*(Be_x*y12 + Be_y*x21));
   	Ru = dueb + BetaGradu + gamma*ueb - feb;
	normRu = fabs(Ru);	
   	Gradu[0] = 0.5*invArea*(y23*ue1 + y31*ue2 + y12*ue3);
   	Gradu[1] = 0.5*invArea*(x32*ue1 + x13*ue2 + x21*ue3);
   	normGradu = sqrt(Gradu[0]*Gradu[0]+Gradu[1]*Gradu[1]);
	normu = fabs(ueb);
	if (normu>=1e-12 && normGradu>=1e-12)
		Eu = 0.25*h*normRu/normGradu + 0.125*h*h*normRu/normu;
	else
		Eu = 0;

	return Eu;
} 



