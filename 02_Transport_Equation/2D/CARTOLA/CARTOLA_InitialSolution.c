#include "../01_CommonFiles/TranspEquation.h"
#include "cartola.h" 

int CARTOLA_InitialSolution(ParametersType *Parameters, NodeType* Node, double *u)
{
	int I, nnodes;
	double X, Y;
	double ro = 0.25, xo = 0.5, yo = 0.5 ;

	nnodes = Parameters->nnodes;

	for(I=0;I<nnodes;I++){
		if (Node[I].Type>0){
			X = Node[I].x;
			Y = Node[I].y;
			if(((X-xo)*(X-xo) - (Y-yo)*(Y-yo)<=ro*ro)&&((X-xo)*(X-xo) - (Y-yo)*(Y-yo)>=0))
				u[Node[I].invP_id] = 0.5;
			else
				u[Node[I].invP_id] = 0;
		}
	}
	
	return 0;
}


