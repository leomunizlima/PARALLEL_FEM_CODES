#include "TranspEquation.h"
#include "../../../00_CommonFiles/Reordering/reordering.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

void insertNode(NodeListType **L, int i, int j, int *nnz_out);
void mount_IA_JA(int n, NodeListType **L,  int *IA, int *JA);

void Permutation_of_LM(ParametersType *Parameters, FemStructsType *FemStructs, int *Perm, int *invPerm)
{
	int nel = Parameters->nel;
	int NEQ = Parameters->NEQ;
	int neq_bef = Parameters->neq_bef;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	

	/********************************************
		Obtain adjacency list of unknowns
	*********************************************/

	int I, J, K, va, vb, nnz = 0;
	NodeListType **List = calloc(NEQ,sizeof(NodeListType*)); 
	
	for (K=0; K<nel; K++){
		for (I=0;I<3;I++){
			va = Element[K].Vertex[I];
			if (Node[va].Type != 0){
				for (J=0; J<3; J++){
					vb = Element[K].Vertex[J];
					if (Node[vb].Type != 0){
						insertNode(List, Node[va].invP_id + neq_bef, Node[vb].invP_id + neq_bef, &nnz);
					}
				}

			}

		}
		
	}

	/***********************************
		Convert to CSR 
	************************************/
	int *IA, *JA, *PermCSR;
	IA = mycalloc("JA of 'Permute_LM'", NEQ+1,sizeof(int));
	JA = mycalloc("JA of 'Permute_LM'", nnz,sizeof(int));  	
	PermCSR = mycalloc("PermCSR of 'csr_Inititalization'",nnz+1,sizeof(int));

	Parameters->nnzero = nnz;

	mount_IA_JA(NEQ, List,IA,JA);

	for (I=0; I<NEQ; I++)
		Perm[I] = I;
	Perm[NEQ] = NEQ;	

	/*********************************************
		Applying reordering to get permutation 
	**********************************************/
	Parameters->bandwidth_bef = MATRIX_bandwidth(NEQ, JA,IA);
	reordering(Parameters, JA, IA, Perm);
	MATRIX_ROW_permutation (NEQ, nnz, JA, IA, Perm, PermCSR);
	MATRIX_COL_permutation (NEQ, nnz, JA, IA, Perm, PermCSR);
	Parameters->bandwidth_aft = MATRIX_bandwidth(NEQ, JA,IA);

	/*******************************************
		Inverse permutation
	*******************************************/
	for (I=0; I<NEQ; I++)
		invPerm[Perm[I]] = I;
	invPerm[NEQ] = NEQ;
	
	free(IA);
	free(JA);
	free(PermCSR);
}

void insertNode(NodeListType **L, int i, int j, int *nnz_out)
{
	if (L[i]==NULL){
		L[i] = calloc(1,sizeof(NodeListType));
		L[i]->J = j;
		L[i]->next = NULL;
		(*nnz_out)++;
		return;	
	}

	if (L[i]->J == j) return;

	NodeListType *current, *new;

	new = calloc(1,sizeof(NodeListType));
	new->J = j;
	new->next = NULL;

	if (L[i]->J > j){
		new->next = L[i];
		L[i] = new;
		(*nnz_out)++;
		return;
	}

	current = L[i];
	while(current->next != NULL && current->next->J < j)
			current = current->next;
	
	if (current->next == NULL){
		current->next = new;
		(*nnz_out)++;
		return;
	}
	
	if (current->next->J == j){
		return;
	}

	new->next = current->next;		
	current->next = new;
	(*nnz_out)++;
}

void mount_IA_JA(int n, NodeListType **L,  int *IA, int *JA)
{
	int I, nnz;
	NodeListType *current, *temp;

	nnz = 0;
	IA[0] = 0;
	for (I=0; I<n; I++){
		current = L[I];
		while (current !=NULL){
			JA[nnz] = current->J;
			temp = current;
			current = current->next;
			free(temp);
			nnz++;
		}	
		IA[I+1] = nnz;
	}
	free(L);

	return;
}


