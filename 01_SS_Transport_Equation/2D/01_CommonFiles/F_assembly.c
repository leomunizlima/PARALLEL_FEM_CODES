#include "SSTranspEquation.h"

void F_assembly(NodeType *Node, int J1, int J2, int J3, double Fe1, double Fe2, double Fe3, double *F, double ke[3][3],double X[3], double Y[3], FemFunctionsType *FemFunctions)
{
	if (Node[J1].Type == 1){ 
		F[Node[J1].invP_id] += Fe1;	
		if (Node[J2].Type == 0)
			F[Node[J1].invP_id] -= ke[0][1]*FemFunctions->upresc(X[1],Y[1]);
		if (Node[J3].Type == 0)
			F[Node[J1].invP_id] -= ke[0][2]*FemFunctions->upresc(X[2],Y[2]);
	}

	if (Node[J2].Type == 1){ 
		F[Node[J2].invP_id] += Fe2;	
		if (Node[J1].Type == 0)
			F[Node[J2].invP_id] -= ke[1][0]*FemFunctions->upresc(X[0],Y[0]);
		if (Node[J3].Type == 0)
			F[Node[J2].invP_id] -= ke[1][2]*FemFunctions->upresc(X[2],Y[2]);
	}

	if (Node[J3].Type == 1){ 
		F[Node[J3].invP_id] += Fe3;	
		if (Node[J1].Type == 0)
			F[Node[J3].invP_id] -= ke[2][0]*FemFunctions->upresc(X[0],Y[0]);
		if (Node[J2].Type == 0)
			F[Node[J3].invP_id] -= ke[2][1]*FemFunctions->upresc(X[1],Y[1]);
	}
}


