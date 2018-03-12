#include "NavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/Reordering/reordering.h"

void csr_Initialization(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	int neq, neq_bef, neq_aft, NEQ, nel, nnodes, nnzero, nnzero_bef, nnzero_aft, count, count2, pos, K, I, J, II, JJ;
	int size = NNOEL*NDOF;
	int size2 = size*size;
	int **CSR_by_Element, **CSR_by_Element_bef, **CSR_by_Element_aft, *IA, *JA;
	NodeListType *current, *temp, **CSR_List, **CSR_List_bef, **CSR_List_aft;
	int **lm = FemStructs->lm;
	NodeType *Node = FemStructs->Node;

	nel  = Parameters->nel;
	neq  = Parameters->neq;
	NEQ = Parameters->NEQ;
	nnodes = Parameters->nnodes;
	neq_bef = Parameters->neq_bef;
	neq_aft = NEQ-neq-neq_bef;

	CSR_List = (NodeListType **)mycalloc("CSR_List of 'CSR_Initialization'",neq,sizeof(NodeListType *));
	for (K = 0; K<neq ; K++)
		CSR_List[K] = NULL;	

	CSR_List_bef = (NodeListType **)mycalloc("CSR_List_bef of 'CSR_Initialization'",neq,sizeof(NodeListType *));
	for (K= 0; K<neq; K++)
		CSR_List_bef[K] = NULL;

	CSR_List_aft = (NodeListType **)mycalloc("CSR_List_aft of 'CSR_Initialization'",neq,sizeof(NodeListType *));
	for (K= 0; K<neq; K++)
		CSR_List_aft[K] = NULL;
	
	nnzero = 0;
	nnzero_bef = 0;
	nnzero_aft = 0;

	for (K = 0; K < nel ; K++)
	{
		for (I=0; I<size; I++){
			II = lm[K][I];
			if (II >= neq_bef && II < neq + neq_bef){
				for (J=0; J<size; J++){
					JJ = lm[K][J];
					if (JJ < neq_bef)
						csr_List_insertA( CSR_List_bef, II-neq_bef,JJ, &nnzero_bef);
					else if (JJ < neq + neq_bef)
						csr_List_insertA( CSR_List, II-neq_bef,JJ-neq_bef, &nnzero);
					else if (JJ != NEQ)
						csr_List_insertA( CSR_List_aft, II-neq_bef,JJ-neq-neq_bef, &nnzero_aft);
				}
			}	


		}

	}

	Parameters->nnzero = nnzero;
	Parameters->nnzero_bef = nnzero_bef;
	Parameters->nnzero_aft = nnzero_aft;

	CSR_by_Element = (int**) mycalloc("CSR_by_Element of 'CSR_Initialization'", nel,sizeof(int*));
	CSR_by_Element_bef = (int**) mycalloc("CSR_by_Element_bef of 'CSR_Initialization'", nel,sizeof(int*));
	CSR_by_Element_aft = (int**) mycalloc("CSR_by_Element_aft of 'CSR_Initialization'", nel,sizeof(int*));

	for (K = 0; K < nel; K++){
		CSR_by_Element[K] = (int*) mycalloc("CSR_by_Element of 'CSR_Initialization'",size2,sizeof(int));
		CSR_by_Element_bef[K] = (int*) mycalloc("CSR_by_Element of 'CSR_Initialization'",size2,sizeof(int));
		CSR_by_Element_aft[K] = (int*) mycalloc("CSR_by_Element of 'CSR_Initialization'",size2,sizeof(int));
	}

	
	for (K = 0; K < nel ; K++){
		for (I=0; I<size2; I++){
			CSR_by_Element[K][I] = nnzero;
			CSR_by_Element_bef[K][I] = nnzero_bef;
			CSR_by_Element_aft[K][I] = nnzero_aft;
		}
	}


	JA = (int *)mycalloc("JA of 'CSR_Initialization'",nnzero, sizeof(int));
	IA = (int *)mycalloc("IA of 'CSR_Initialization'",neq+1, sizeof(int));
	count = 0;
	for (K = 0 ; K < neq ; K++){
		current = CSR_List[K];
		while (current != NULL){
			JA[count++] = current->J;	
			IA[K+1]++;
			current = current->next;
		}	
	}	 

	for (K = 1; K < neq; K++)
		IA[K+1] += IA[K];	
	
	
	for (K = 0; K < nel ; K++){
		
		for (I=0; I<size; I++){
	
			II = lm[K][I];

			if (II >= neq_bef && II < neq + neq_bef){

				for (J=0; J<size; J++){

					JJ = lm[K][J];

					if (JJ >= neq_bef && JJ < neq + neq_bef){

						pos = csr_search( II-neq_bef, JJ-neq_bef, CSR_List[II-neq_bef]);

						CSR_by_Element[K][I*size+J] = IA[II-neq_bef] + pos;
					}		
				}

			}
			
	
		}

	}

	for (K = 0; K < neq; K++){
		current = CSR_List[K];
		while (current != NULL){
			temp = current;		
			current	= current->next;
			free(temp);
		}
	}

	free(CSR_List);

	
	if (neq_bef>0){
		
		int neqaux_bef;
		
		count2 = 0;
		for (K = 0; K < neq; K++)
			if (CSR_List_bef[K] != NULL)
				count2++;

		neqaux_bef = count2; // Count the number of lines in Struct before has at least one nonzero coefficient
		Parameters->neqaux_bef = neqaux_bef;		

		int *JA_bef = (int *)mycalloc("JA_bef of 'CSR_Initialization'",nnzero_bef, sizeof(int));
		int *IA_bef = (int *)mycalloc("IA_bef of 'CSR_Initialization'",neqaux_bef+1, sizeof(int));
		int *aux_bef = (int *)mycalloc("aux_bef of 'CSR_Initialization'",neqaux_bef, sizeof(int));
		int *IAidex_bef = (int *)mycalloc("IAide_bef of 'CSR_Initialization",neq+1,sizeof(int));

		MatrixData->JA_bef = JA_bef;			
		MatrixData->IA_bef = IA_bef;			
		MatrixData->aux_bef = aux_bef; 

		count = 0;
		count2 = 0;

		for (K = 0 ; K <= neq; K++)
			IAidex_bef[K] = neqaux_bef;

		for (K = 0 ; K < neq ; K++){
			current = CSR_List_bef[K];
			if (current != NULL){
				IAidex_bef[K] = count2;
				aux_bef[count2] = K;
				count2++;
			}
			while (current != NULL){
				JA_bef[count++] = current->J;	
				IA_bef[count2]++;
				current = current->next;
			}
		}	
		for (K = 1; K < neqaux_bef; K++)
			IA_bef[K+1] += IA_bef[K];	

		for (K = 0; K < nel ; K++){
			
			for (I=0; I<size; I++){

				II = lm[K][I];

				if (II >= neq_bef && II < neq + neq_bef){
		
					for (J=0; J<size; J++){

						JJ = lm[K][J];
				
						if (JJ < neq_bef){
							pos = csr_search( II-neq_bef, JJ, CSR_List_bef[II-neq_bef]);
							CSR_by_Element_bef[K][I*size+J] = IA_bef[IAidex_bef[II-neq_bef]] + pos;
						}
					}
				}
		
			}
	
		}
		
		free(IAidex_bef);

	}

	for (K = 0; K < neq; K++){
		current = CSR_List_bef[K];
		while (current != NULL){
			temp = current;		
			current	= current->next;
			free(temp);
		}
	}
	
	free(CSR_List_bef);


	if (neq_aft>0){
	
		int neqaux_aft;
	
		count2 = 0;
		for (K = 0; K < neq; K++)
			if (CSR_List_aft[K] != NULL)
				count2++;

		neqaux_aft = count2; // Count the number of lines in Struct after has at least one nonzero coefficient
		Parameters->neqaux_aft = neqaux_aft;		

		int *JA_aft = (int *)mycalloc("JA_bef of 'CSR_Initialization'",nnzero_aft, sizeof(int));
		int *IA_aft = (int *)mycalloc("IA_bef of 'CSR_Initialization'",neqaux_aft+1, sizeof(int));
		int *aux_aft = (int *)mycalloc("IA_bef of 'CSR_Initialization'",neqaux_aft, sizeof(int));
		int *IAidex_aft = (int *)mycalloc("IAide_bef of 'CSR_Initialization",neq+1,sizeof(int));

		MatrixData->JA_aft = JA_aft;			
		MatrixData->IA_aft = IA_aft;			
		MatrixData->aux_aft = aux_aft; 

		count = 0;
		count2 = 0;

		for (K = 0 ; K <= neq; K++)
			IAidex_aft[K] = neqaux_aft;

		for (K = 0 ; K < neq ; K++){
			current = CSR_List_aft[K];
			if (current != NULL){
				aux_aft[count2] = K;
				IAidex_aft[K] = count2;
				count2++;
			}
			while (current != NULL){
				JA_aft[count++] = current->J;	
				IA_aft[count2]++;
				current = current->next;
			}	
		}
		for (K = 1; K < neqaux_aft; K++)
			IA_aft[K+1] += IA_aft[K];	

		for (K = 0; K < nel ; K++){

			for (I=0 ; I<size; I++){
	
				II = lm[K][I];
	
				if (II >= neq_bef && II < neq + neq_bef){
			
					for (J=0; J<size; J++){	
						
						JJ = lm[K][J];
	
						if (JJ>=neq+neq_bef && JJ != NEQ){
							pos = csr_search( II-neq_bef, JJ-neq-neq_bef, CSR_List_aft[II-neq_bef]);
							CSR_by_Element_aft[K][I*size+J] = IA_aft[IAidex_aft[II-neq_bef]] + pos;
						}
					}
		
				}
	
			}

		}

		free(IAidex_aft);

	}
	
	for (K = 0; K < neq; K++){
		current = CSR_List_aft[K];
		while (current != NULL){
			temp = current;		
			current	= current->next;
			free(temp);
		}
	}
	
	free(CSR_List_aft);
	

	int *PermCSR = mycalloc("PermCSR of 'csr_Inititalization'",nnzero,sizeof(int));
	int *perm = mycalloc("perm of 'csr_Initialization'",neq,sizeof(int));

	for (I=0; I<nnzero; I++)
		PermCSR[I] = I;
	for (I=0; I<neq; I++)
		perm[I]=I;

	Parameters->bandwidth_bef = MATRIX_bandwidth(neq, JA,IA);
	reordering (Parameters,JA,IA,perm);
	MATRIX_ROW_permutation (neq, nnzero, JA, IA, perm, PermCSR);
	MATRIX_COL_permutation (neq, nnzero, JA, IA, perm, PermCSR);
	Parameters->bandwidth_aft = MATRIX_bandwidth(neq, JA,IA);

	int *invPermCSR = mycalloc("invPermCSR of 'csr_Inititalization'",nnzero+1,sizeof(int));
	int *invperm = mycalloc("invPermCSR of 'csr_Inititalization'",neq+1,sizeof(int));

	for (I=0;I<nnzero; I++)
		invPermCSR[PermCSR[I]] = I;
	invPermCSR[nnzero]=nnzero;

	for (I=0;I<neq;I++)
		invperm[perm[I]]=I;
	invperm[neq]=neq;

	int aux[size2];
	for (K=0; K<nel; K++){
		for (I=0;I<size2;I++)
			aux[I] = invPermCSR[CSR_by_Element[K][I]];
		for (I=0;I<size2;I++)
			CSR_by_Element[K][I] = aux[I];
		for (I=0; I<size; I++)
			if (lm[K][I]>=neq_bef && lm[K][I]<neq_bef+neq)
				lm[K][I] = invperm[lm[K][I]-neq_bef]+neq_bef;
	}
	
	for (K=0; K<nnodes; K++){
		for (I=0; I<NDOF; I++){
			if (Node[K].Type[I]==1)	
				Node[K].invP_id[I] = invperm[Node[K].invP_id[I]];
		}
	}

	free(PermCSR);
	free(invPermCSR);

	int *tempvet;

	if (neq_bef>0){
			
		int neqaux_bef = Parameters->neqaux_bef;		
		int nsend_bef = Parameters->nsend_bef;		
		int *aux_bef = MatrixData->aux_bef; 
		int *JA_bef = MatrixData->JA_bef;
		int *IA_bef = MatrixData->IA_bef;
		int *IdSend_bef = FemStructs->IdSend_bef;
			
		PermCSR = mycalloc("invPermCSR of 'csr_Inititalization'",nnzero_bef, sizeof(int));
		invPermCSR = mycalloc("invPermCSR of 'csr_Inititalization'",nnzero_bef+1, sizeof(int));
		tempvet = mycalloc("tempvet of 'csr_Inititalization'",neqaux_bef, sizeof(int));

		for (I=0; I<nnzero_bef; I++)
			PermCSR[I] = I;
		
		for (I = 0; I<neqaux_bef; I++)
			tempvet[I] = invperm[aux_bef[I]];
		for (I = 0; I<neqaux_bef; I++)
			aux_bef[I] = tempvet[I];
	
		free(tempvet);

		MATRIX_SPARSE_ROW_permutation (neqaux_bef, nnzero_bef, JA_bef, IA_bef, aux_bef, PermCSR);

		tempvet = mycalloc("tempvet of 'csr_Inititalization'",nsend_bef, sizeof(int));
		
		for (I=0; I<nsend_bef; I++)
			tempvet[I] = invperm[IdSend_bef[I]];
		for (I=0; I<nsend_bef; I++)
			IdSend_bef[I] = tempvet[I];
			
		free(tempvet);

		for (I=0;I<nnzero_bef; I++)
			invPermCSR[PermCSR[I]] = I;
		invPermCSR[nnzero_bef] = nnzero_bef;

		for (K=0; K<nel; K++){
			for (I=0;I<size2;I++)
				aux[I] = invPermCSR[CSR_by_Element_bef[K][I]];
			for (I=0;I<size2;I++)
				CSR_by_Element_bef[K][I] = aux[I];
		}

	
		free(PermCSR);
		free(invPermCSR);

	}
	if (neq_aft>0){

		int neqaux_aft = Parameters->neqaux_aft;		
		int nsend_aft = Parameters->nsend_aft;		
		int *aux_aft = MatrixData->aux_aft; 
		int *JA_aft = MatrixData->JA_aft;
		int *IA_aft = MatrixData->IA_aft;
		int *IdSend_aft = FemStructs->IdSend_aft;

		PermCSR = mycalloc("invPermCSR of 'csr_Inititalization'",nnzero_aft, sizeof(int));
		invPermCSR = mycalloc("invPermCSR of 'csr_Inititalization'",nnzero_aft+1, sizeof(int));
		tempvet = mycalloc("tempvet of 'csr_Inititalization'",neqaux_aft, sizeof(int));

		for (I=0; I<nnzero_aft; I++)
			PermCSR[I] = I;
		
		for (I = 0; I<neqaux_aft; I++)
			tempvet[I] = invperm[aux_aft[I]];
		for (I = 0; I<neqaux_aft; I++)
			aux_aft[I] = tempvet[I];
	
		free(tempvet);

		MATRIX_SPARSE_ROW_permutation (neqaux_aft, nnzero_aft, JA_aft, IA_aft, aux_aft, PermCSR);

		tempvet = mycalloc("tempvet of 'csr_Inititalization'",nsend_aft, sizeof(int));
		
		for (I=0; I<nsend_aft; I++)
			tempvet[I] = invperm[IdSend_aft[I]];
		for (I=0; I<nsend_aft; I++)
			IdSend_aft[I] = tempvet[I];
			
		free(tempvet);

		for (I=0; I<nnzero_aft; I++)
			invPermCSR[PermCSR[I]] = I;
		invPermCSR[nnzero_aft] = nnzero_aft;
	
		for (K=0; K<nel; K++){
			for (I=0;I<size2;I++)
				aux[I] = invPermCSR[CSR_by_Element_aft[K][I]];
			for (I=0;I<size2;I++)
				CSR_by_Element_aft[K][I] = aux[I];
		}

		free(PermCSR);
		free(invPermCSR);
	}

	MatrixData->Perm = perm;	
	MatrixData->Scheme_by_Element = CSR_by_Element;
	MatrixData->Scheme_by_Element_bef = CSR_by_Element_bef;
	MatrixData->Scheme_by_Element_aft = CSR_by_Element_aft;
	MatrixData->JA = JA; 
	MatrixData->IA = IA;
	free(invperm);	
}


void csr_List_insertA(NodeListType **CSR_List, int I, int J, int *nnzero_out)
{
	int nnzero;
	NodeListType *current, *previous, *new;

	nnzero = *nnzero_out;
	nnzero++;
	new = mycalloc("new of 'CSR_Initialization'", 1, sizeof(NodeListType));
	new->J = J;
	new->next = NULL;
	//only to initialization
	previous = new;
	previous->next = NULL;

	if (CSR_List[I]==NULL)
		CSR_List[I] = new;
	else if (CSR_List[I]->J > J){
		new->next = CSR_List[I];
		CSR_List[I] = new;
	}
	else {
		current = CSR_List[I];		
		while (current != NULL){
			if (J <= current->J)			
				break;
			previous = current;	
			current = current->next;
		}
		if (current == NULL)
			previous->next = new;
		else if (current->J == J){
			nnzero--;
			free(new);
		}
		else {
			new->next = current;
			previous->next = new;
		}
	}
	*nnzero_out = nnzero;
}

int csr_search(int I, int J, NodeListType *current)
{
	int position;

	position = -1;

	while (current != NULL){
		position++;
		if (current->J == J)
			break;	 
		current = current->next;

	}

	return position;
}



