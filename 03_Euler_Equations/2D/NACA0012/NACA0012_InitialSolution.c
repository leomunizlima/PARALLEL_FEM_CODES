# include "naca0012.h"

int NACA0012_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u){

	int I, nnodes;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++)
	{
		if (Node[I].Type[0] != 0){ u[Node[I].invP_id[0]] = 1.0; }
		if (Node[I].Type[1] != 0){ u[Node[I].invP_id[1]] = 1.0; }
		if (Node[I].Type[2] != 0){ u[Node[I].invP_id[2]] = 0.0; }
		if (Node[I].Type[3] != 0){ u[Node[I].invP_id[3]] = 10.5; }
	}

	return 0;
}


