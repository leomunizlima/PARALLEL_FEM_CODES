#include "EulerEquations.h"

void dimensionless_Ax_Ay_calculations(double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4])
{
	// *** Jacobian matrices Ax and Ay : primitive variables
	double v1 = Ub[1] / Ub[0];
	double v2 = Ub[2] / Ub[0];
	double E  = Ub[3] / Ub[0];
	double norma_v = v1*v1 + v2*v2; //||v||^2

	/*Dimensionless Ax and Ay*/
	double epsilon, epsilon2;
	epsilon = Mach; //reference Mach number
	epsilon2 = epsilon * epsilon;

	// *** Ax coefficients
	Ax[0][0] = 0.0;
	Ax[0][1] = 1.0;
	Ax[0][2] = 0.0;
	Ax[0][3] = 0.0;
	Ax[1][0] = -(v1 * v1) + ((gamma - 1) * norma_v * 0.5) ;
	Ax[1][1] = (3.0 - gamma) * v1;
	Ax[1][2] = -(gamma - 1.0) * v2;
	Ax[1][3] = (gamma - 1.0)/epsilon2;
	Ax[2][0] = - v1*v2;
	Ax[2][1] = v2;
	Ax[2][2] = v1;
	Ax[2][3] = 0.0;
	Ax[3][0] = - v1 * (gamma * E - ((gamma - 1.0) * epsilon2 * norma_v));
	Ax[3][1] = gamma * E - (gamma - 1.0) * epsilon2 * ( 0.5 * norma_v + v1 * v1 );
	Ax[3][2] = -(gamma - 1.0) * epsilon2 * v1 * v2;
	Ax[3][3] = gamma * v1;

	// *** Ay coefficients
	Ay[0][0] = 0.0;
	Ay[0][1] = 0.0;
	Ay[0][2] = 1.0;
	Ay[0][0] = 0.0;
	Ay[1][0] = Ax[2][0];
	Ay[1][1] = v2;
	Ay[1][2] = v1;
	Ay[1][3] = 0.0;
	Ay[2][0] = -(v2 * v2) + (((gamma - 1.0) * norma_v * 0.5)); 
	Ay[2][1] = -(gamma - 1.0) * v1;
	Ay[2][2] = (3.0 - gamma) * v2;
	Ay[2][3] = Ax[1][3];
	Ay[3][0] = - v2 * (gamma * E - (gamma - 1.0) * epsilon2 * norma_v);
	Ay[3][1] = Ax[3][2];
	Ay[3][2] = gamma * E - (gamma - 1.0) * epsilon2 * (0.5 * norma_v + v2 * v2);
	Ay[3][3] = gamma * v2;

}

void dimensional_Ax_Ay_calculations(double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4])
{
	double rho = Ub[0];
	double v1 = Ub[1] / Ub[0];
	double v2 = Ub[2] / Ub[0];
	double E  = Ub[3] / Ub[0];
	double norma_U23 = (Ub[1] * Ub[1]) + (Ub[2] * Ub[2]);

	// *** Ax coefficients
	Ax[0][0] = 0.0;
	Ax[0][1] = 1.0;
	Ax[0][2] = 0.0;
	Ax[0][3] = 0.0;
	Ax[1][0] = -(v1 * v1) + (((gamma - 1) * norma_U23 * 0.5) / (rho * rho));
	Ax[1][1] = (3.0 - gamma) * v1;
	Ax[1][2] = -(gamma - 1.0) * v2;
	Ax[1][3] = gamma - 1.0;
	Ax[2][0] = - v1*v2;
	Ax[2][1] = v2;
	Ax[2][2] = v1;
	Ax[2][3] = 0.0;
	Ax[3][0] = (-gamma * v1 * E) + ((gamma - 1.0) * v1 * (norma_U23 / (rho * rho)));
	Ax[3][1] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v1 * v1)));
	Ax[3][2] = -((gamma - 1.0) * v1 * v2);
	Ax[3][3] = gamma * v1;

	// *** Ay coefficients

	Ay[0][0] = 0.0;
	Ay[0][1] = 0.0;
	Ay[0][2] = 1.0;
	Ay[0][0] = 0.0;
	Ay[1][0] = Ax[2][0];
	Ay[1][1] = v2;
	Ay[1][2] = v1;
	Ay[1][3] = 0.0;
	Ay[2][0] = -(v2 * v2) + (((gamma - 1.0) * norma_U23 * 0.5) / (rho * rho)); 
	Ay[2][1] = -(gamma - 1.0) * v1;
	Ay[2][2] = (3.0 - gamma) * v2;
	Ay[2][3] = Ax[1][3];
	Ay[3][0] = (- gamma * v2 * E) + (((gamma - 1.0) * v2 * norma_U23) / (rho * rho));
	Ay[3][1] = Ax[3][2];
	Ay[3][2] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v2 * v2)));
	Ay[3][3] = gamma * v2;
}	

