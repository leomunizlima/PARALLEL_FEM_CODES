#include "EulerEquations.h"

int uB_InitialSolution(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *u, double *uB)
{
	int i, e, nel, NEQ, eNDOF, J1, J2, J3;
	double x1, x2, x3, y1, y2, y3;
	double third = 1.0/3.0;
	double Ue[12];
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	nel = Parameters->nel;
	NEQ = Parameters->NEQ;

	for (e = 0; e < nel; e++){
		// *** Global node that composes the element
		J1 = Element[e].Vertex[0];
		J2 = Element[e].Vertex[1];
		J3 = Element[e].Vertex[2];
	
		// *** Nodal coordinates and differential operator
		x1 = Node[J1].x;
		x2 = Node[J2].x;
		x3 = Node[J3].x;
		
		y1 = Node[J1].y;
		y2 = Node[J2].y;
		y3 = Node[J3].y;

		// *** Fill Ue with the values of the previous solution or workaround for each degree of liberty of the element node
		// Density
	 	if (Node[J1].invP_id[0] != NEQ)
	   		Ue[0] = u[Node[J1].invP_id[0]];
	   	else
	   		Ue[0] = FemFunctions->rhopresc(x1,y1);
	   	
	   	if (Node[J2].invP_id[0] != NEQ)
			Ue[4] = u[Node[J2].invP_id[0]];
		else
			Ue[4] = FemFunctions->rhopresc(x2,y2);
	   
		if (Node[J3].invP_id[0] != NEQ)
			Ue[8] = u[Node[J3].invP_id[0]];
		else
			Ue[8] = FemFunctions->rhopresc(x3,y3);
		
		// Velocity x
		if (Node[J1].invP_id[1] != NEQ)
	   		Ue[1] = u[Node[J1].invP_id[1]];
	   	else
	   		Ue[1] = FemFunctions->v1presc(x1, y1);
	   
		if (Node[J2].invP_id[1] != NEQ)
			Ue[5] = u[Node[J2].invP_id[1]];
		else
			Ue[5] = FemFunctions->v1presc(x2, y2);
		
	   	if (Node[J3].invP_id[1] != NEQ)
			Ue[9] = u[Node[J3].invP_id[1]];
		else
			Ue[9] = FemFunctions->v1presc(x3, y3);
	
		// Velocity y
		if (Node[J1].invP_id[2] != NEQ)
	   		Ue[2] = u[Node[J1].invP_id[2]];
	   	else
	   		Ue[2] = FemFunctions->v2presc(x1, y1);
	   
		if (Node[J2].invP_id[2] != NEQ)
			Ue[6] = u[Node[J2].invP_id[2]];
		else
			Ue[6] = FemFunctions->v2presc(x2, y2);
	   
		if (Node[J3].invP_id[2] != NEQ)
			Ue[10] = u[Node[J3].invP_id[2]];
		else
			Ue[10] = FemFunctions->v2presc(x3, y3);
	
		// Energy
		if (Node[J1].invP_id[3] != NEQ)
	   		Ue[3] = u[Node[J1].invP_id[3]];
	   	else
	   		Ue[3] = FemFunctions->epresc(x1, y1);
	   	
	   	if (Node[J2].invP_id[3] != NEQ)
			Ue[7] = u[Node[J2].invP_id[3]];
		else
			Ue[7] = FemFunctions->epresc(x2, y2);
	   
		if (Node[J3].invP_id[3] != NEQ)
			Ue[11] = u[Node[J3].invP_id[3]];
		else
			Ue[11] = FemFunctions->epresc(x3, y3);


		// *** The triangle centroid
		eNDOF = e*NDOF;
		for (i = 0; i < 4; i++){
			uB[eNDOF+i] = (Ue[i] + Ue[i+4] + Ue[i+8]) * third;
		}

	}	
	return 0;
}// end build


