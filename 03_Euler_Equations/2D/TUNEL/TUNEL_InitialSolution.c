#include "tunel.h"

int TUNEL_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u)
{
	int I, nnodes;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++)
	{
		if (Node[I].Type[0] != 0){ // U1 = rho incognita
			u[Node[I].invP_id[0]] = 1.4;
		}
		if (Node[I].Type[1] != 0){ // U2 = rho*v1 = 1.4 * 3.0 incognita
			u[Node[I].invP_id[1]] = 4.2;
		}
		if (Node[I].Type[2] != 0){ // U3 = rho*v2 = 1.4 * 0.0 incognita
			u[Node[I].invP_id[2]] = 0.0;
		}
		if (Node[I].Type[3] != 0){ // U4 = rho*e = 1.4 * 6.285714286 incognita
			u[Node[I].invP_id[3]] = 8.8;
		}
	}

	return 0;
}

