#include "SSNavierStokesEquations.h"

int Fill_LM(ParametersType *Parameters, int **lm, NodeType *Node, ElementType *Element)
{
	int e, i, j, k, l;
	int NEQ = Parameters->NEQ;
	int nel = Parameters->nel;
	int neq_bef = Parameters->neq_bef;


	for (e = 0; e < nel; e++)   //'for' what through all the elements
	{
		l = 0;   //'l' assists in the assembly lm
		for (i = 0; i < NNOEL; i++)   //'for' what through all the nodes of elements
		{
			j = Element[e].Vertex[i];   //'j' stores the global value of the node
			for(k = 0; k < NDOF; k++)   //'for' what through all the degrees of freedom
			{
				if(Node[j].Type[k] == 0){
					lm[e][l] = NEQ;
					l++;
				}else{
					lm[e][l] = Node[j].invP_id[k] + neq_bef;
					l++;
				}
			}
		}
	}

	return 0;
}


