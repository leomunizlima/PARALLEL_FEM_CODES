# include "obliquo.h"

int OBLIQUO_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u){
	int I, nnodes;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++)
	{
		if (Node[I].Type[0] != 0){ u[Node[I].invP_id[0]] = 1.0; }
		if (Node[I].Type[1] != 0){ u[Node[I].invP_id[1]] =  0.98481; }//0.984807753; }
		if (Node[I].Type[2] != 0){ u[Node[I].invP_id[2]] = -0.17365; }//-0.173648177; }
		if (Node[I].Type[3] != 0){ u[Node[I].invP_id[3]] = 0.94643; }// 0.946425; }
	}

	return 0;
}


