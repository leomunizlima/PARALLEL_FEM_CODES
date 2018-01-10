#include "TranspEquation.h"

void eval_U_dU(ParametersType *Parameters,FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *U,double *dU)
{
	int I;
	double x, y;
	int nnodes = Parameters->nnodes;
	double *u = FemStructs->u;
	double *du = FemStructs->du;
	NodeType *Node = FemStructs->Node;

	for (I=0; I<nnodes; I++){

		x = Node[I].x;
		y = Node[I].y;

		if (Node[I].Type == 0){
			U[I] = FemFunctions->upresc(x,y);
			dU[I] = 0.0;
		}
		else{
			U[I] = u[Node[I].invP_id];
			dU[I] = du[Node[I].invP_id];
		}
/*		if (Node[I].Type == 1){
			U[I] = u[Node[I].invP_id];
			dU[I] = du[Node[I].invP_id];
		}
		else if (Node[I].Type > 1){
			U[I] = 0;//0.333*u[Node[I].invP_id];
			dU[I] = 0;//0.333*du[Node[I].invP_id];
		}
		else{
			U[I] = FemFunctions->upresc(x,y);
			dU[I] = 0.0;
		}*/
	}
}



