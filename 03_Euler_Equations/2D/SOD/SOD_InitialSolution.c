#include "sod.h"

int SOD_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u)
{
	int I, nnodes;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++)
	{
		if(Node[I].x < 0.5){ // se o nó está a esquerda da membrana
			if (Node[I].Type[0] != 0){ // U1 = rho incognita
				u[Node[I].invP_id[0]] = 1.0;
			}
			if (Node[I].Type[1] != 0){ // U2 = rho*v1 = 1.0 * 0.0 incognita
				u[Node[I].invP_id[1]] = 0.0;
			}
			if (Node[I].Type[2] != 0){ // U3 = rho*v2 = 1.0 * 0.0 incognita
				u[Node[I].invP_id[2]] = 0.0;
			}
			if (Node[I].Type[3] != 0){ // U4 = rho*e = 1.0 * 2.5 incognita
				u[Node[I].invP_id[3]] = 2.5;
			}
		}else if(Node[I].x >= 0.5){ // se o nó está a direita da membrana
			if (Node[I].Type[0] != 0){ // U1 = rho incognita
				u[Node[I].invP_id[0]] = 0.125;
			}
			if (Node[I].Type[1] != 0){ // U2 = rho*v1 = 0.125 * 0.0 incognita
				u[Node[I].invP_id[1]] = 0.0;
			}
			if (Node[I].Type[2] != 0){ // U3 = rho*v2 = 0.125 * 0.0 incognita
				u[Node[I].invP_id[2]] = 0.0;
			}
			if (Node[I].Type[3] != 0){ // U4 = rho*e = 0.125 * 2.0 incognita
				u[Node[I].invP_id[3]] = 0.25;
			}
		}
	}

	return 0;
}

