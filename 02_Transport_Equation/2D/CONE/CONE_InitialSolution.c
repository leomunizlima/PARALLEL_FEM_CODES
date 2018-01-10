#include "cone.h"

int CONE_InitialSolution(ParametersType *Parameters, NodeType* Node, double *u)
{
	int I, nnodes;
	double X, Y;
	double R;

	nnodes = Parameters->nnodes;

	for(I=0;I<nnodes;I++){
		if (Node[I].Type!=0){
			X = Node[I].x;
			Y = Node[I].y;
			R = (X-5)*(X-5)+ (Y-7.5)*(Y-7.5);
			u[Node[I].invP_id] = exp(-0.5*R);
		}
	}

	return 0;
}


