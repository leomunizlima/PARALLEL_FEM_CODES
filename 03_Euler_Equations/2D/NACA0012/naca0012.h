# include "../01_CommonFiles/EulerEquations.h"

double NACA0012_rhopresc(double, double);
double NACA0012_v1presc(double, double);
double NACA0012_v2presc(double, double);
double NACA0012_epresc(double, double);
double NACA0012_gamma(double, double);
double NACA0012_cv(double, double);
double NACA0012_theta(double, double);
int NACA0012_InitialSolution(ParametersType *, NodeType *, double *);
void NACA0012_BC_no_penetrability(int, int, int, NodeType *, double, double, double, double, double, double, double, double, double, double, 
				  double [4][4], double [4][4], double [12][12], double [12][4], double [4][12], double [12], double [4], double [12], double [12], double [4], double [4]);



