#include"protos.h"
#include <time.h>

#define NNOEL 3

int main(int argc, char **argv)
{
	int neq, neq_mesh, nnodes, nel, I, J, K, II, JJ, NDOF, tag;
	char FileName[300];
	FILE *InFile;
	NodeType *Node;
	ElementType *Element;
	struct timespec Start, End;
	double Elapsed_Time;
	
	

	/**************************************************************************************************/
	//					Testing initial parameters
	/**************************************************************************************************/
	if (argc<6)
	{
		printf("Use ./Preprocess <Mesh.dat> <Problem> <NDOF> <Reordering> <path to create RANK files> <list of number of subdivisions>\n");
		exit(1);
	}


	/*****************************************************************************/
	//				Reading nodes	
	/*****************************************************************************/		
	sprintf(FileName,"%s",argv[1]);
	InFile = fopen(FileName, "r");
	tag = fscanf(InFile, "%d", &nnodes);
	NDOF = atoi(argv[3]);
	Node = (NodeType*) calloc(nnodes, sizeof(NodeType));
	int *Node_mesh = calloc(nnodes, sizeof(int));
	for (I=0, neq = 0, neq_mesh =0; I<nnodes; I++)
	{
		tag = fscanf(InFile, "%lf%lf", &(Node[I].x), &(Node[I].y));

		for (J=0; J<NDOF; J++){
			tag = fscanf(InFile,"%d", &(Node[I].Type[J]));
			if (Node[I].Type[J] == 1)
				Node[I].id[J] = neq++;
			else
				Node[I].id[J] = -1;
		}
		if (find_neq_mesh(Node[I].id,NDOF))
			Node_mesh[I] = neq_mesh++;
		else
			Node_mesh[I] = -1;

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

	if (tag<0){
		printf ("Error in some parameter\n");
		exit(1);
	}
	/*************************************************************************************************************/
	
	/********************************************
		Obtain adjacency list of unknowns
	*********************************************/

	int va, vb, nnz = 0, nnz_mesh = 0;
	ListType **List = calloc(neq,sizeof(ListType*)); 
	ListType **List_mesh = calloc(neq_mesh,sizeof(ListType*)); 
	
	clock_gettime(CLOCK_MONOTONIC, &Start);
	for (K=0; K<nel; K++){
		for (I=0;I<3;I++){
			va = Element[K].Vertex[I];
			for (II = 0 ; II<NDOF; II++){
				if (Node[va].Type[II] != 0){
					for (J=0; J<3; J++){
						vb = Element[K].Vertex[J];
						for (JJ=0; JJ<NDOF; JJ++){
							if (Node[vb].Type[JJ] != 0)
								insertNode(List, Node[va].id[II], Node[vb].id[JJ], &nnz);
						}
					}

				}
			}
			if (Node_mesh[va] != -1){
				for (J=0; J<3; J++){
					vb = Element[K].Vertex[J];
					if (Node_mesh[vb] != -1){
						insertNode(List_mesh, Node_mesh[va], Node_mesh[vb], &nnz_mesh);
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
	int *IA, *JA, *IA_mesh, *JA_mesh;
	IA = calloc(neq+1,sizeof(int));
	JA = calloc(nnz,sizeof(int));  	
	mount_IA_JA(neq,nnz,List,IA,JA);
	IA_mesh = calloc(neq_mesh+1,sizeof(int));
	JA_mesh = calloc(nnz_mesh,sizeof(int));  	
	mount_IA_JA(neq_mesh,nnz_mesh,List_mesh,IA_mesh,JA_mesh);

	/*********************************
		Applying reordering
	*********************************/

	clock_gettime(CLOCK_MONOTONIC, &Start);
	int *P = calloc(neq,sizeof(int));
	int *P_mesh = calloc(neq_mesh,sizeof(int));
	if (strcasecmp(argv[4],"Spectral")==0) //Spectral
		REORDERING_SPECTRAL (neq_mesh, nnz_mesh, JA_mesh,  IA_mesh, P_mesh);
	else if (strcasecmp(argv[4],"SYMRCM")==0) //Symtric_RCM
		REORDERING_SYMRCM (neq_mesh, nnz_mesh,JA_mesh,  IA_mesh, P_mesh);
	else if (strcasecmp(argv[4],"NOT")==0){ //No reordering
		for (I=0; I<neq_mesh; I++)
			P_mesh[I] = I;	
	}
	else{
		printf("Reordering scheme is not defined correctly!\n");
		exit(1);
	}
	free(JA_mesh);
	free(IA_mesh);
	
	int a;
	int **P_aux = calloc(neq_mesh,sizeof(int*));
	int *invP_mesh = calloc(neq_mesh,sizeof(int*));
	
	for (I=0; I<neq_mesh; I++){
		P_aux[I] = calloc(NDOF,sizeof(int));
		invP_mesh[P_mesh[I]]=I;
	}
	free(P_mesh);
		
	for (I=0; I<nnodes; I++){
		if (Node_mesh[I] != -1){
			a = Node_mesh[I];
			for (J=0; J<NDOF; J++)
				P_aux[invP_mesh[a]][J] = Node[I].id[J];	
		}
	}
	free(invP_mesh);
	free(Node_mesh);

	for (I=0,neq=0; I<neq_mesh; I++){
		for (J=0;J<NDOF;J++){
			if (P_aux[I][J] != -1)
				P[neq++] = P_aux[I][J];
		}
		free(P_aux[I]);
	}
	free(P_aux);

	printf("Bandwidth before reordering = %d\n", MATRIX_bandwidth(neq,IA, JA));

	MATRIX_ROW_permutation (neq, nnz, IA, JA, P);
	MATRIX_COL_permutation (neq, nnz, IA, JA, P);
	
	printf("Bandwidth after reordering = %d\n", MATRIX_bandwidth(neq,IA, JA));

	clock_gettime(CLOCK_MONOTONIC, &End);
	
	Elapsed_Time = End.tv_sec - Start.tv_sec + 1e-9*(End.tv_nsec - Start.tv_nsec);
	printf("Runtime of %s reordering was %lf seconds.\n", argv[4], Elapsed_Time);

	int *invP = calloc(neq,sizeof(int));
		for (I=0; I<neq; I++)
			invP[P[I]] = I;


//	int *M=calloc(neq*neq,sizeof(int));
/*	for (I=0; I<neq; I++)
		for (J = IA[I]; J<IA[I+1]; J++)
			M[I*neq+JA[J]] = 1;
	for (I=0; I<neq; I++){
		for (J=0; J<neq; J++)
			printf("%d ",M[I*neq + J]);
		printf("\n");
	}*/
	/***********************************
		Partitioning
	**********************************/
	int N, nP, *d, *L;
	int number_of_executions = argc - 6;

	for (N=0; N<number_of_executions; N++){
		
		nP = atoi(argv[N+6]);	

		PARTITIONING_MIN_MAX (neq, IA, JA, nP, &d, &L);

		//*********************************************************************

		/**********************************
			Dividing elements
		**********************************/
		ListType **ElemByPart; 
		int v, count;
		
		ElemByPart = calloc(nel,sizeof(ListType*));
		count = 0;
		for (I=0; I<nel; I++){
			for (J=0; J<3; J++){
				v = Element[I].Vertex[J];
				for (JJ=0; JJ<NDOF; JJ++){
					if (Node[v].Type[JJ] != 0){
						for (K=0; K<nP; K++){
							if(invP[Node[v].id[JJ]]>=d[K] && invP[Node[v].id[JJ]]<d[K+1]){
								insertNode(ElemByPart, I, K, &count);
	
							}
						}			
					}
				}
			}	
		}

		
			
		/**********************************
			Dividing nodes
		**********************************/	
		ListType **NodeByPart;
	
		NodeByPart = calloc(nnodes,sizeof(ListType *));

		count = 0;
		for (I=0; I<nel; I++){
			for (K=0; K<nP; K++){
				if(Search(ElemByPart,I,K)){
					for (J=0; J<3; J++){
						v = Element[I].Vertex[J];
						insertNode(NodeByPart, v, K, &count);
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
	
		Problem = argv[2];
		map = calloc(nnodes, sizeof(int));
	
		for (K=0; K<nP; K++){
			sprintf(filename,"RANK_%02d_%s_%d_%d.dat",K,Problem,nnodes,nel);
			f = fopen(filename,"w");
			count_nodes = COUNTER(NodeByPart,nnodes, K);
			fprintf(f,"%d\n",count_nodes);
			count_nodes = 0;
			for (I=0; I<nnodes; I++){

				if (Search(NodeByPart,I,K)){
					fprintf(f,"%.14lf\t%.14lf\t", Node[I].x, Node[I].y);
					for (JJ=0; JJ<NDOF; JJ++){
						Type = catch_Type(Node,I,K,JJ, d,invP);
						if (K>0 && Search(NodeByPart,I,K-1)) 
							Type_bef = catch_Type(Node,I,K-1,JJ,d,invP);
						else
							Type_bef = 1;
						if (K<nP-1 && Search(NodeByPart,I,K+1))
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
				if (Search(ElemByPart,I,K)){
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
			printf("RANK %d: %.2lf%% (Percentual of overlapping)\n", K, 100*(nP*count_elem/(double)nel-1.0)); 
			
			
		}

		char command[300];

		command[0]='\0';

		sprintf(command,"mkdir -p %s/%02d\n",argv[5],nP); 
		tag = system(command);
		if (tag != 0){
			printf("Executation error in command 1!\n");
			exit(1);
		} 
		
		command[0]='\0';
		sprintf(command,"mv RANK*.dat %s/%02d\n",argv[5],nP); 
		tag = system(command);
		if (tag != 0){
			printf("Executation error in command 2!\n");
			exit(1);
		}
		
		for (I=0; I<nP+1; I++)
			printf("d[%d]=%d\n",I,d[I]);

		free(map);
		free(d);
		free(L);
		free_List(ElemByPart,nel);
		free_List(NodeByPart,nnodes);
	}

	free(Node);
 	free(Element);		
	free(P);
	free(invP);
	free(IA);
	free(JA);
	
	return 0;
}

int COUNTER(ListType **L, int n, int K)
{
	int I, count = 0;

	for(I=0;I<n;I++){
		if (Search(L,I,K))
			count++;
	}
	return count;
}

int Search(ListType **L, int I, int K)
{
	ListType *Temp = L[I];

	while(Temp!=NULL && Temp->J !=K)
		Temp = Temp->next;
	
	if (Temp != NULL)
		return 1;
	else
		return 0;		
}

void free_List(ListType **L, int n)
{
	int I;
	ListType *current, *temp;

	for (I=0; I<n; I++){
		current = L[I];
		while (current !=NULL){
			temp = current;
			current = current->next;
			free(temp);
		}	
	}
	free(L);

	return;

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
		free(new);
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

int find_neq_mesh(int *vec, int n)
{
	int i=0;
	int found=0;

	for (i=0;i<n;i++){
		if (vec[i]!=-1){
			found = 1;
			break;
		}
	}

	return found;
}

int MATRIX_bandwidth (int n, int *IA, int *JA)
{
	int i;
	int bandl, bandr;
	int bandwidth=0;
	
	bandl = 0;
	bandr = 0;
	for (i = 0; i < n ; i++)
	{
		if (fabs(i - JA[IA[i]]) > bandl)
			bandl= i - JA[IA[i]];
		if (fabs(JA[IA[i+1]-1]-i) > bandr);	
			bandr = JA[IA[i+1]-1]-i; 
	}
	if (bandl>bandr)
		bandwidth = bandl;
	else
		bandwidth = bandr;

	return bandwidth;
}





