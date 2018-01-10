# include "../01_CommonFiles/EulerEquations.h"

double CYLINDER_rhopresc(double, double);
double CYLINDER_v1presc(double, double);
double CYLINDER_v2presc(double, double);
double CYLINDER_epresc(double, double);
double CYLINDER_gamma(double, double);
double CYLINDER_cv(double, double);
double CYLINDER_theta(double,double);
int CYLINDER_InitialSolution(ParametersType *, NodeType *, double *);
void CYLINDER_BC_no_penetrability(int, int, int, NodeType *, double, double, double, double, double, double, double, double, double, double, 
				  double [4][4], double [4][4], double [12][12], double [12][4], double [4][12], double [12], double [4], double [12], double [12], double [4], double [4]);


