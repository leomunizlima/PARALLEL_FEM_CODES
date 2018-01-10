#include "baroclinic.h"

int BAROCLINIC_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u)
{
	int I, nnodes;
	double rho0 = 1.0, rho1 = 0.001, rho2 = 1.8, p0 = 1.0, gamma = 1.4;
	double rho, p1, p, Phiy, x, y, u0, epsilon, epsilon2, v1, U2, U4;
	nnodes = Parameters->nnodes;
	epsilon = Parameters->Mach;
	u0 = sqrt(gamma);
	p1 = gamma;
	epsilon2 = epsilon * epsilon;
	
	for(I = 0; I < nnodes; I++)
	{
		x = Node[I].x;
		y = Node[I].y;	

		//function rho
		if( y <= 4.0 )
			Phiy = (rho2/8.0) * y;
		else
			Phiy = rho2 * (y/8.0 - 1.0);	
		rho = rho0 + 0.5 * rho1 * epsilon * ( 1.0 + cos( PI * x / 20.0 ) ) + Phiy;	
		//function v1
		v1 = 0.5 * u0 * ( 1.0 + cos( PI * x / 20.0 ) );
		//function v2 = 0.0
		
		//function E
		p = p0 + 0.5 * epsilon * p1 * ( 1.0 + cos( PI * x / 20.0 ) );
		U2 = v1 * rho;		
		U4 = p / ( gamma - 1.0 ) + 0.5 * epsilon2 * rho * ( v1 * v1 ); //state equation and U4 = rho*E. ||v||^2 = v1*v1, v2=0.
		

		if (Node[I].Type[0] > 0){ // U1 = rho incognita
			u[Node[I].invP_id[0]] = rho;
		}
		if (Node[I].Type[1] > 0){ // U2 = rho*v1 
			u[Node[I].invP_id[1]] = U2;
		}
		if (Node[I].Type[2] > 0){ // U3 = rho*v2 
			u[Node[I].invP_id[2]] = 0.0;
		}
		if (Node[I].Type[3] > 0){ // U4 = rho*E
			u[Node[I].invP_id[3]] = U4;
		}
	}

	return 0;
}

