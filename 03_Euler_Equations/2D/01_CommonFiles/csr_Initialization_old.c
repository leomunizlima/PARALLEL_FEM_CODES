#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

void List_insert(NodeListType **List, int I, int J, int *counter);
void Fill_Id_NodeData(NodeListType *List, NodeDataType *NodeData);

void csr_Initialization(int nnodes, int nel, int *nnz_out, NodeType *Node, ElementType *Element, NodeListType **List, CSRAuxType *CSRAux)
{
	int I, J, I1, I2, I3, nnz = 0;
	int *counter = (int*) mycalloc("counter of 'csr_Initialization'",nnodes,sizeof(int));
	NodeDataType **NodeData = (NodeDataType**) mycalloc("NodeData of 'csr_Inititalization'",nnodes,sizeof(NodeDataType*));

	for (I=0; I<nel; I++)
	{
 		I1 = Element[I].Vertex[0];
		I2 = Element[I].Vertex[1];
		I3 = Element[I].Vertex[2];
 
		List_insert(List, I1, I1, counter);
		List_insert(List, I1, I2, counter);
		List_insert(List, I1, I3, counter);
		List_insert(List, I2, II, counter);
		List_insert(List, I2, I2, counter);
		List_insert(List, I2, I3, counter);
		List_insert(List, I3, II, counter);
		List_insert(List, I3, I2, counter);
		List_insert(List, I3, I3, counter);
	}

	for (I=0; I<nnodes; I++){
		NodeData[I] = (NodeDataType*) mycalloc("NodeData of 'csr_Initialization'",counter[I],sizeof(NodeDataType));
		Fill_Id_NodeData(List[I],NodeData[I]);	
	}

	for (I= 0; I < nnodes; I++){
		current = List[I];
		while (current != NULL){
			temp = current;		
			current	= current->next;
			free(temp);
		}

	}
	free(List);

	for (I= 0; I<nnodes; I++){
		for (J=0; J<counter[I]; J++){
			for (K=0; K<NDOF; K++){
				if (Node[NodeData[I][J].Id].id[K]>=0)
					nnz++;
			}
		}
	}

	CSRAux->conter = counter;
	CSRAux->NodeData = NodeData;
	
}


void List_insert(NodeListType **List, int I, int J, int *counter)
{
	NodeListType *current=NULL, *previous=NULL, *new=NULL;

	new = mycalloc("new of 'ede_Initialization'", 1,sizeof(NodeListType));
	new->J = J;
	new->next = NULL;
	counter[I]++;

	if (List[I]==NULL)
		List[I] = new;
	else if (List[I]->J > J){
		new->next = List[I];
		List[I] = new;
	}
	else {
		current = List[I];		
		while (current != NULL){
			if (J <= current->J)			
				break;
			previous = current;	
			current = current->next;
		}
		if (current == NULL)
			previous->next = new;
		else if (current->J == J){
			free(new);
			counter[I]--;
		}
		else {
			new->next = current;
			previous->next = new;
		}
	}
}

void Fill_Id_NodeData(NodeListType *List, NodeDataType *NodeData)
{
	J = 0;

	while(List!=NULL){
		NodeData[J++] = List->J;
		List = List->next;
	}

}

