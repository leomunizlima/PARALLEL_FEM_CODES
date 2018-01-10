#include "Mesh_Reordering.h"
/**********************************************************************************************************/
//					Reordering using triangular RCM					
/**********************************************************************************************************/

void reordering_rcm(int nnos, int nel, int NDOF, char *Problem, double *X, double *Y, int **Type, ElementType *Element)
{
	char command[300], filename[300];
	FILE *In,*Out;
	int I, J;
	int *perm = mycalloc("perm",nnos,sizeof(int));
	int *perm_inv = mycalloc("perm_inv",nnos,sizeof(int));
	int *tempInt = mycalloc("temp",nnos,sizeof(int));
	double *tempDoub = mycalloc("temp",nnos,sizeof(double));
		

	strcpy(filename,"");
	sprintf(filename,"%s_nodes.txt",Problem);
	Out = myfopen(filename,"w");
	for (I=0;I<nnos; I++)
		fprintf(Out,"%lf %lf\n",X[I],Y[I]); 
	fclose(Out);

	strcpy(filename,"");
	sprintf(filename,"%s_elements.txt",Problem);
	Out = myfopen(filename,"w");
	for (I=0; I<nel; I++){
		fprintf(Out, "%d %d %d\n", Element[I].Vertex[0],  Element[I].Vertex[1], Element[I].Vertex[2]);
	}	
	fclose(Out);
	
	strcpy(command,"");
	sprintf(command, "echo '%s' > %s", Problem,Problem);
	if (!system(command)) printf("Ops!\n");

	strcpy(command,"");
	sprintf(command, "./triangulation_rcm < %s", Problem);
	if (!system(command)) printf("\nFinalizing rcm triangulation...\n");
	
	strcpy(filename,"");
	sprintf(filename,"%s_rcm_perm.txt", Problem);
	In = myfopen(filename,"r");
	
	for (I=0; I<nnos; I++){
		if (!fscanf(In,"%d", &perm[I])) printf("Ops!\n");
		perm[I]--;
	}
	fclose(In);


		
	for (I=0; I<nnos; I++)
		perm_inv[perm[I]] = I;

	for (I=0; I<nnos; I++)
		tempDoub[I] = X[perm[I]];
	
	for (I=0; I<nnos; I++)
		X[I] = tempDoub[I];

	for (I=0; I<nnos; I++)
		tempDoub[I] = Y[perm[I]];
	
	for (I=0; I<nnos; I++)
		Y[I] = tempDoub[I];

	for (J=0 ; J<NDOF; J++){
		
		for (I=0; I<nnos; I++)
			tempInt[I] = Type[perm[I]][J];
	
		for (I=0; I<nnos; I++)
			Type[I][J] = tempInt[I];

	}

	int Node;	
	for (I=0; I<nel; I++){
		Node = Element[I].Vertex[0];
		Element[I].Vertex[0] = perm_inv[Node];
		
		Node = Element[I].Vertex[1];
		Element[I].Vertex[1] = perm_inv[Node];
		
		Node = Element[I].Vertex[2];
		Element[I].Vertex[2] = perm_inv[Node];
	}		
	//Removing auxiliar files 
	strcpy(command,"\0");
	sprintf(command, "rm -rf %s*.txt",Problem);
	if (!system(command))printf("removing temporary files...\n");	

	return;
}

