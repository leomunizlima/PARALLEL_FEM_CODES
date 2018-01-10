#include"protos.h"
#include <time.h>

#define NNOEL 3

int main(int argc, char **argv)
{
	int neq, nnodes, nnzero, nel, I, J, K, JJ, NDOF, tag, size;
	double *F, *u, *Diag;
	char FileName[300], label[300], Reordering[300];
	FILE *InFile;
	NodeType *Node;
	ElementType *Element;
	struct timespec Start, End;
	double Elapsed_Time;
	
	

	/**************************************************************************************************/
	//					Testing initial parameters
	/**************************************************************************************************/
	if (argc!=6)
	{
		printf("Use ./Preprocess <Mesh.dat> <nP> <Problem> <NDOF> <Reordering>\n");
		exit(1);
	}


	/*****************************************************************************/
	//				Reading nodes	
	/*****************************************************************************/		
	sprintf(FileName,"%s",argv[1]);
	InFile = fopen(FileName, "r");
	tag = fscanf(InFile, "%d", &nnodes);
	NDOF = atoi(argv[4]);
	Node = (NodeType*) calloc(nnodes, sizeof(NodeType));
	for (I=0, neq = 0; I<nnodes; I++)
	{
		tag = fscanf(InFile, "%lf%lf", &(Node[I].x), &(Node[I].y));

		for (J=0; J<NDOF; J++){
			tag = fscanf(InFile,"%d", &(Node[I].Type[J]));
			if (Node[I].Type[J] == 1)
				Node[I].id[J] = neq++;
			else
				Node[I].id[J] = -1;
		}
	}
	/*****************************************************************************/

	
	/************************************************************************************************************/	
	//						Reading connection mesh
	/************************************************************************************************************/	
	tag = fscanf(InFile, "%d", &nel);
	Element = (ElementType*) calloc(nel, sizeof(ElementType));
	for (I=0; I<nel; I++){
		tag = fscanf(InFile, "%d%d%d", &(Element[I].Vertex[0]),  &(Element[I].Vertex[1]), &(Element[I].Vertex[2]));
		if (NDOF>1)
			tag = fscanf(InFile, "%d", &(Element[I].Type));
	}
	fclose(InFile);
	/*************************************************************************************************************/
	
	/********************************************
		Obtain adjacency list of unknowns
	*********************************************/

	int va, vb, nnz = 0;
	int nP = atoi(argv[2]);	
	ListType **List = calloc(neq,sizeof(ListType*)); 
	
	clock_gettime(CLOCK_MONOTONIC, &Start);
	for (K=0; K<nel; K++){
		for (I=0;I<3;I++){
			va = Element[K].Vertex[I];
			for (JJ = 0 ; JJ<NDOF; JJ++){
				if (Node[va].Type[JJ] != 0){
					for (J=0; J<3; J++){
						vb = Element[K].Vertex[J];
						if (Node[vb].Type[JJ] != 0){
							insertNode(List, Node[va].id[JJ], Node[vb].id[JJ], &nnz);
						}
					}

				}
			}

		}
		
	}
	clock_gettime(CLOCK_MONOTONIC, &End);
	
	Elapsed_Time = End.tv_sec - Start.tv_sec + 1e-9*(End.tv_nsec - Start.tv_nsec);
	printf("Runtime to create adjacency list was %lf seconds.\n", Elapsed_Time);


	/***********************************
		Convert to CSR 
	************************************/
	int *IA, *JA;
	IA = calloc(neq+1,sizeof(int));
	JA = calloc(nnz,sizeof(int));  	
	mount_IA_JA(neq,nnz,List,IA,JA);

	/*********************************
		Applying reordering
	*********************************/

	clock_gettime(CLOCK_MONOTONIC, &Start);
	
	int *P = calloc(neq,sizeof(int));
	if (strcasecmp(argv[5],"Spectral")==0) //Spectral
		REORDERING_SPECTRAL (neq, nnz, JA,  IA, P);
	else if (strcasecmp(argv[5],"SYMRCM")==0) //Symtric_RCM
		REORDERING_SYMRCM (neq, nnz,JA,  IA, P);
	else{
		printf("Reordering scheme is not defined correctly!\n");
		exit(1);
	}

	MATRIX_ROW_permutation (neq, nnz, IA, JA, P);
	MATRIX_COL_permutation (neq, nnz, IA, JA, P);
	
	clock_gettime(CLOCK_MONOTONIC, &End);
	
	Elapsed_Time = End.tv_sec - Start.tv_sec + 1e-9*(End.tv_nsec - Start.tv_nsec);
	printf("Runtime of %s reordering was %lf seconds.\n", argv[5], Elapsed_Time);
	

	/***********************************
		Partitioning
	**********************************/
	int *d, *L;
	PARTITIONING_MIN_MAX (neq, IA, JA, nP, &d, &L);


	/**********************************
		Dividing elements
	**********************************/
	int **ElemByPart, *ElemByPartAux; 
	int v;

	ElemByPart = calloc(nel,sizeof(int*));
	ElemByPartAux = calloc(nel*nP,sizeof(int)); 
	for (I=0; I<nel; I++)
		ElemByPart[I] = &ElemByPartAux[I*nP];
	
	for (I=0; I<neq; I++)
		invP[P[I]] = I;
	
	for (I=0; I<nel; I++){
		for (J=0; J<3; J++){
			v = Element[I].Vertex[J];
			for (JJ=0; JJ<NDOF; JJ++){
				if (Node[v].Type[JJ] != 0){
					for (K=0; K<nP; K++){
						if(invP[Node[v].id[JJ]]>=d[K] && invP[Node[v].id[JJ]]<d[K+1]){
							ElemByPart[I][K] = 1;
						}
					}			
				}
			}
		}	
	}
			
	/**********************************
		Dividing nodes
	**********************************/	
	int **NodeByPart, *NodeByPartAux;
	
	NodeByPart = calloc(nnodes,sizeof(int*));
	NodeByPartAux = calloc(nnodes*nP,sizeof(int));
	for (I=0; I<nnodes; I++)
		NodeByPart[I] = &NodeByPartAux[I*nP];
	
	for (I=0; I<nel; I++){
		for (J=0; J<3; J++){
			v = Element[I].Vertex[J];
			for (K=0; K<nP; K++){
				if(ElemByPart[I][K]==1){
					NodeByPart[v][K] = 1;
				}			
			}
		}	
	}
	
	/********************************
		Dividing files
	*********************************/
	int count_nodes, count_elem, Type, Type_bef, Type_aft, *map;
	char filename[100], *Problem;	
	FILE *f;
	
	Problem = argv[3];
	map = calloc(nnodes, sizeof(int));
	
	for (K=0; K<nP; K++){
		sprintf(filename,"RANK_%02d_%s_%d_%d.dat",K,Problem,nnodes,nel);
		f = fopen(filename,"w");
		count_nodes = COUNTER(NodeByPart,nnodes, K);
		fprintf(f,"%d\n",count_nodes);
		count_nodes = 0;
		for (I=0; I<nnodes; I++){
			if (NodeByPart[I][K] == 1){
				fprintf(f,"%.14lf\t%.14lf\t", Node[I].x, Node[I].y);
				for (JJ=0; JJ<NDOF; JJ++){
					Type = catch_Type(Node,I,K,JJ, d,invP);
					if (K>0 && NodeByPart[I][K-1]==1) 
						Type_bef = catch_Type(Node,I,K-1,JJ,d,invP);
					else
						Type_bef = 1;
					if (K<nP-1 && NodeByPart[I][K+1]==1)
						Type_aft = catch_Type(Node, I,K+1,JJ,d,invP);
					else
						Type_aft = 1;
					if (Type == 0)
						fprintf(f,"%d\t", Type); 
					else{
						if (Type==1 && Type_bef==3 && Type_aft==2)
							fprintf(f,"%d\t23\t%d\t",Type,invP[Node[I].id[JJ]]-d[K]); 
						else if (Type==1 && Type_bef==3)
							fprintf(f,"%d\t2\t%d\t",Type,invP[Node[I].id[JJ]]-d[K]); 
						else if (Type==1 && Type_aft==2)
							fprintf(f,"%d\t3\t%d\t",Type,invP[Node[I].id[JJ]]-d[K]); 
						else
							fprintf(f,"%d\t1\t%d\t",Type, invP[Node[I].id[JJ]]-d[K]); 
							
					}
				}
				fprintf(f,"\n");
				map[I] = count_nodes++;
			}
		}	 
		count_elem = COUNTER(ElemByPart,nel,K);
		fprintf(f,"%d\n",count_elem);
		for (I=0; I<nel; I++){
			if (ElemByPart[I][K]==1){
				for (J=0; J<3; J++){
					v = Element[I].Vertex[J];
					fprintf(f,"%d ", map[v]);
				}
				if (NDOF>1)
					fprintf(f,"%d\n",Element[I].Type);
				else
					fprintf(f,"\n");
			}
		}
		fclose(f);
	}

	for (I=0; I<nP+1; I++)
		printf("d[%d]=%d\n",I,d[I]);

	for (I=0;I<neq;I++)
		printf("%d\n",P[I]);


	free(map);
	free(Node);
 	free(Element);		
	free(ElemByPartAux);
	free(ElemByPart);
	free(NodeByPartAux);
	free(NodeByPart);
	free(d);
	free(L);
	free(P);
	free(invP);
	free(IA);
	free(JA);
	
	return 0;
}

int COUNTER(int **M, int n, int K)
{
	int I, count = 0;

	for(I=0;I<n;I++){
		if (M[I][K]==1)
			count++;
	}
	return count;
}

int catch_Type(NodeType *Node, int I, int K, int JJ, int *d, int *invP)
{
	int Type;
	if (Node[I].Type[JJ]==0)
		Type = 0;
	else if (invP[Node[I].id[JJ]]<d[K])
		Type = 2;
	else if (invP[Node[I].id[JJ]]>=d[K+1])
		Type = 3;
	else
		Type = 1;

	return Type;
}
void insertNode(ListType **L, int i, int j, int *nnz_out)
{
	if (L[i]==NULL){
		L[i] = calloc(1,sizeof(ListType));
		L[i]->J = j;
		L[i]->next = NULL;
		(*nnz_out)++;
		return;	
	}

	if (L[i]->J == j) return;

	ListType *current, *new;

	new = calloc(1,sizeof(ListType));
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

void mount_IA_JA(int n, int nnz, ListType **L,  int *IA, int *JA)
{
	int I;
	ListType *current, *temp;

	
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

/*----------------------------------------------------------------------------
 * Perform the row permutation
 *--------------------------------------------------------------------------*/
void MATRIX_ROW_permutation (int n, int nnz, int *IA, int *JA, int* p)
{
	int i, j, k;

	int*    auxJA = calloc( nnz  ,sizeof(int));
	int*    auxIA = calloc((n+1) ,sizeof(int));
  
	auxIA[0] = 0;
	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = IA[p[i]]; j <= IA[p[i]+1] - 1; ++j)
		{
			auxJA[k] = JA[j];
			k  = k + 1;
		}
		auxIA[i+1] = k;    
	}

	memcpy(&JA[0],&auxJA[0],nnz*sizeof(int));
	memcpy(&IA[0],&auxIA[0],(n+1)*sizeof(int));

	free(auxJA);
	free(auxIA);
}



/*----------------------------------------------------------------------------
* Perform the colunm permutation
*--------------------------------------------------------------------------*/
void MATRIX_COL_permutation (int n, int nnz, int *IA, int *JA, int* p)
{
	int i, j, k;

	ARRAY* a = calloc (nnz,sizeof(ARRAY));
	int*   q = calloc (n ,sizeof(int));
	
	for (i = 0; i < n; ++i) 
		q[p[i]] = i; 

	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = IA[i]; j <= IA[i+1] - 1; ++j)
		{
			a[k].arr2 = q[JA[j]];
			a[k].arr3 = i;
				k = k + 1;
		}
			IA[i+1] = k;    
	}

	qsort(a,nnz,sizeof(ARRAY),COMPARE_array);
	
	for (i = 0; i < nnz; ++i)
		JA[i] = a[i].arr2;
	

	free(a);
	free(q);
}

int COMPARE_array (const void * a, const void * b)
{
	if (((ARRAY*)a)->arr3 <  ((ARRAY*)b)->arr3) return -1;
	if (((ARRAY*)a)->arr3 >  ((ARRAY*)b)->arr3) return  1;
	if (((ARRAY*)a)->arr3 == ((ARRAY*)b)->arr3)
	{
		if (((ARRAY*)a)->arr2 < ((ARRAY*)b)->arr2) return -1;
		if (((ARRAY*)a)->arr2 > ((ARRAY*)b)->arr2) return  1;
	}
	return 0;
}





