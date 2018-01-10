#include "TranspEquation.h"

double CAU_DD_ShockCapture(double kx, double ky, double Be_x, double Be_y, double gamma, double ue1, double ue2, double ue3, double ueb, double dueb, double feb, 
                        double y23, double y31, double y12, double x32, double x13, double x21, double invArea, double h)
{
	double cbtil, BetaGradu, Ru, Gradu[2], normGradu;	
	double tolEps = 1e-10;

	BetaGradu = 0.5*invArea*(ue1*(Be_x*y23 + Be_y*x32) + ue2*(Be_x*y31 + Be_y*x13) + ue3*(Be_x*y12 + Be_y*x21));
   	Ru = dueb + BetaGradu + gamma*ueb - feb;
   	Gradu[0] = 0.5*invArea*(y23*ue1 + y31*ue2 + y12*ue3);
   	Gradu[1] = 0.5*invArea*(x32*ue1 + x13*ue2 + x21*ue3);
   	normGradu = sqrt(Gradu[0]*Gradu[0]+Gradu[1]*Gradu[1]);
 	if (normGradu > tolEps)
		cbtil = 0.5*h*fabs(Ru)/normGradu;
	else
		cbtil = 0;

	return cbtil;
}


