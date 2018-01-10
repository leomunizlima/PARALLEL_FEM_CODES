#include "../01_CommonFiles/TranspEquation.h"
#include "pudim.h" 

int PUDIM_InitialSolution(ParametersType *Parameters, NodeType* Node, double *u)
{
	int I, nnodes;
	double X, Y;
	
	nnodes = Parameters->nnodes;

	for(I=0; I<nnodes; I++){
		if (Node[I].Type>0){
			X = Node[I].x;
			Y = Node[I].y;
			if (fabs(X)<=1e-14 && Y<=0.)
				u[Node[I].invP_id] = -sin(PI*Y);
			else
				u[Node[I].invP_id] = 0;
		}
	}
	return 0;
}


