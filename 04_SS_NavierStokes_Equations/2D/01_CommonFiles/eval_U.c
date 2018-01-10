#include "SSNavierStokesEquations.h"

void eval_U(ParametersType *Parameters,FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *U)
{
	int I;
	double x, y;
	int nnodes = Parameters->nnodes;
	double *u = FemStructs->u;
	NodeType *Node = FemStructs->Node;

	for (I=0; I<nnodes; I++){
		
		x = Node[I].x;
		y = Node[I].y;
	
		//velocity u
		if (Node[I].Type[0]>0)
	   		U[3*I] = u[Node[I].invP_id[0]];
	   	else
	   		U[3*I]= FemFunctions->v1presc(x,y);
	   	
		//velocity v
		if (Node[I].Type[1]>0)
	   		U[3*I+1] = u[Node[I].invP_id[1]];
	   	else
	   		U[3*I+1] = FemFunctions->v2presc(x,y);

	   	//pression p
		if (Node[I].Type[2]>0)
	   		U[3*I+2]= u[Node[I].invP_id[2]];
	   	else
	   		U[3*I+2] = FemFunctions->ppresc(x,y);
	   	
	}
}
