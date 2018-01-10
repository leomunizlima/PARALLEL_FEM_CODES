#include "SSTranspEquation.h"

int Fill_LM(ParametersType *Parameters, int **lm, NodeType *Node, ElementType *Element)
{
	int J,e,JJ;
	int NEQ = Parameters->NEQ;
	int nel = Parameters->nel;
	int neq_bef = Parameters->neq_bef;

	for (e=0; e<nel; e++)
	{
		lm[e][0]=NEQ;
		lm[e][1]=NEQ;
		lm[e][2]=NEQ;
		for (J=0; J<3; J++)
		{
			JJ = Element[e].Vertex[J];

			if (Node[JJ].Type != 0)
				lm[e][J] = Node[JJ].invP_id + neq_bef;

		}
	}

	return 0;

}


