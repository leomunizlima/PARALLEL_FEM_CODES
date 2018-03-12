#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct
{
	double x,y;
}NodeType;

typedef struct
{
	int Vertex[3];
}ElementType;

int print_bound[16][4]= {{1,1,1,1}/*mark 0 */, {1,1,1,0}/*mark 1 */, {1,1,0,1}/*mark 2 */, {1,1,0,0}/*mark 3 */,
				{1,0,1,1}/*mark 4 */, {1,0,1,0}/*mark 5 */, {1,0,0,1}/*mark 6 */, {1,0,0,0}/*mark 7 */,
				{0,1,1,1}/*mark 8 */, {0,1,1,0}/*mark 9 */, {0,1,0,1}/*mark 10*/, {0,1,0,0}/*mark 11*/,
				{0,0,1,1}/*mark 12*/, {0,0,1,0}/*mark 13*/, {0,0,0,1}/*mark 14*/, {0,0,0,0}/*mark 15*/};

int setmark_according_boundary(int nnodes, double **coord, int *mark, char *problem);

int main(int argc, char **argv)
{
	int nel, nnodes, I, index, aux1, aux2, aux3, aux4, aux5, tag=0;
	double z;	
	NodeType *Node;
	ElementType *Element;
	char label[300];
	char *aux;
	const char delimiter[3]=" \n";
	FILE *In;


	if (argc<3){
		printf("Use ./gmsh2Mesh <msh file> <problem>\n");
		exit(1);
	}

	if ((In=fopen(argv[1],"r"))==NULL){
		printf("File %s not found!\n", argv[1]);
		exit(1);
	}		
	
	//Reading header information
	do{
		tag = fscanf(In,"%s",label);
	}while(strcmp(label, "$Nodes")!=0);

	tag = fscanf(In,"%d",&nnodes);
	
	if ((Node = calloc(nnodes,sizeof(NodeType)))==NULL){
		printf("Memmory allocation error in Node!\n");
		exit(1);
	}

	//Reading nodes
	for (I=0; I<nnodes; I++)
		tag = fscanf(In,"%d%lf%lf%lf", &index, &Node[I].x, &Node[I].y, &z);
	
	//Reading header information
	do{
		tag = fscanf(In,"%s",label);
	}while(strcmp(label, "$Elements")!=0);

	tag = fscanf(In,"%d\n",&nel);

	char *error;
	do{
		error = fgets(label, 300, In);
		aux = strtok(label, delimiter);
		aux = strtok(NULL, delimiter);
		nel--;
	}while (strcmp(aux, "2")!=0);

	nel++;
	
	if ((Element = calloc(nel,sizeof(ElementType)))==NULL){
		printf("Memmory allocation error in Element!\n");
		exit(1);
	}

	aux = strtok(NULL, delimiter);
	aux = strtok(NULL, delimiter);
	aux = strtok(NULL, delimiter);
	aux = strtok(NULL, delimiter);
	Element[0].Vertex[0] = atoi(aux);
	aux = strtok(NULL, delimiter);
	Element[0].Vertex[1] = atoi(aux);
	aux = strtok(NULL, delimiter);
	Element[0].Vertex[2] = atoi(aux);
	
	for(I=1; I<nel;I++)
		tag = fscanf(In, "%d%d%d%d%d%d%d%d", &aux1, &aux2, &aux3, &aux4, &aux5, &Element[I].Vertex[0], &Element[I].Vertex[1], &Element[I].Vertex[2]);  	
	
	fclose(In);
	
	if (tag<0 || error==NULL){
		printf ("Error in some parameter\n");
		exit(1);
	}
	


	//setting mark according to boundary conditions
	double **coord;
	int i, **ien, *mark;
	char *problem;

	coord = calloc(nnodes,sizeof(double*));
	ien = calloc(nel,sizeof(int*));
	mark = calloc(nnodes,sizeof(int));
	problem = argv[2];	


	for (i=0; i<nnodes; i++)
		coord[i] = calloc(2,sizeof(double));
	
	for (i=0; i<nel; i++)
		ien[i] = calloc(3,sizeof(int));

	for (i=0; i<nnodes; i++){
		coord[i][0] = Node[i].x;
		coord[i][1] = Node[i].y;
	}	
	free(Node);

	for (i=0; i<nel; i++){
		ien[i][0] = Element[i].Vertex[0]-1;
		ien[i][1] = Element[i].Vertex[1]-1;
		ien[i][2] = Element[i].Vertex[2]-1;

	}
	free(Element);

	int j,NDOF;
	char FileName[300];
	FILE *File;
	NDOF = setmark_according_boundary(nnodes, coord, mark, problem);

	//Output file 
	sprintf(FileName,"%s_%d_%d.dat", problem, nnodes, nel);
	File = fopen(FileName, "w");
	if (File == NULL)
	{
		printf("File %s not found!\n", FileName);
		exit(1);
	}

	fprintf(File, "%d\n", nnodes);
	for(i = 0; i < nnodes; i++){
		fprintf(File, "%.14lf\t%.14lf\t", coord[i][0], coord[i][1]);
		for (j=0; j<NDOF; j++)
			fprintf(File,"%d\t",print_bound[mark[i]][j]);
		fprintf(File,"\n");
	}

	fprintf(File, "%d\n", nel);
	if (NDOF == 1){
		for(i = 0; i < nel; i++)
			fprintf(File, "%d\t%d\t%d\n", ien[i][0], ien[i][1], ien[i][2]);
	}
	else{
		for(i = 0; i < nel; i++)
			fprintf(File, "%d\t%d\t%d\t-1\n", ien[i][0], ien[i][1], ien[i][2]);
	}
	

	fclose(File);

	//Desalocations
	for (i=0; i<nnodes; i++)
		free(coord[i]);
	for (i=0; i<nel; i++)
		free(ien[i]);
	free(coord);
	free(ien);	
	free(mark);


	return 0;
}

int setmark_according_boundary(int nnodes, double **coord, int *mark, char *problem)
{
	int i, NDOF;

	if(strcasecmp(problem,"OBLIQUO") == 0){
	
		NDOF = 4;
		for (i=0;i<nnodes;i++){
			if (fabs(coord[i][1])<=1e-15)
				mark[i] = 2;
			else if (fabs(coord[i][0])<=1e-15 || fabs(coord[i][1]-1.0)<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}

	}else if(strcasecmp(problem,"REFLETIDO") == 0){
	
		NDOF = 4;
		for (i=0;i<nnodes;i++){
			if (fabs(coord[i][1])<=1e-15)
				mark[i] = 2;
			else if (fabs(coord[i][0])<=1e-15 || fabs(coord[i][1]-1.0)<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}

	}else if(strcasecmp(problem,"EXPLOSION") == 0){

		NDOF = 4;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0])<=1e-15 || fabs(coord[i][1])<=1e-15 || fabs(coord[i][0]-2.0)<=1e-15 || fabs(coord[i][1]-2.0)<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}

	}else if(strcasecmp(problem,"BAROCLINIC") == 0){
	
		NDOF = 4;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][1])<=1e-15 || fabs(coord[i][1]-8.0)<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}

	}else if(strcasecmp(problem,"TUNEL")==0){
		
		NDOF = 4;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][1])<=1e-15 || (fabs(coord[i][1]-0.2)<=1e-15 && coord[i][0]>=0.6) || fabs(coord[i][1]-1.0)<=1e-15)
				mark[i] = 2;
			else if (fabs(coord[i][0]-0.6)<=1e-15 && fabs(coord[i][1]-0.2)<=1e-4)
				mark[i] = 6;
			else if (fabs(coord[i][0]-0.6)<=1e-15 && coord[i][1]<=0.2)
				mark[i] = 4;
			else if (fabs(coord[i][0])<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}

	}else if(strcasecmp(problem,"CYLINDER")==0){
		
		print_bound[6][1]--;
		print_bound[6][2]--;
		NDOF = 4;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]*coord[i][0] + coord[i][1]*coord[i][1] - 0.25)<=1e-10){
				mark[i] = 6;
			}else if (fabs(coord[i][0]+3.0)<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}
	
	}else if(strcasecmp(problem,"NACA0012")==0){
		
		NDOF = 4;
		double x, y, z;
		print_bound[6][1]--;
		print_bound[6][2]--;
 
		for(i=0;i<nnodes;i++){

			x = coord[i][0];
			if (x>0){
				y = 0.178140*sqrt(x) - 0.075600*x - 0.210960*x*x +  0.170580*x*x*x - 0.060900*x*x*x*x;  
				z = -y;
			}else{
				y = 1e5;
				z = 1e5;
			}
			

			if (fabs(x + 0.5)<=1e-8)
				mark[i] = 15;
			else if (fabs(coord[i][1]-y)<=1e-5 || fabs(coord[i][1]-z)<=1e-5){
				mark[i] = 6;
			}else
				mark[i] = 0;
		}

	}else if(strcasecmp(problem,"PUDIM") == 0){ 

		NDOF = 1;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]-1.0)<=1e-15 || fabs(coord[i][0]+1.0)<=1e-15 || fabs(coord[i][1]-1.0)<=1e-15 ||fabs(coord[i][1]+1.0)<=1e-15 ||
				(fabs(coord[i][0])<=1e-15 && fabs(coord[i][1]<=0)))
				mark[i] = 15;
			else
				mark[i] = 0;
		}
	}else if(strcasecmp(problem,"BIGPUDIM") == 0){ 

		NDOF = 1;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]-10.0)<=1e-15 || fabs(coord[i][0]+10.0)<=1e-15 || fabs(coord[i][1]-10.0)<=1e-15 ||fabs(coord[i][1]+10.0)<=1e-15 ||
				(fabs(coord[i][0])<=1e-15 && fabs(coord[i][1]<=0)))
				mark[i] = 15;
			else
				mark[i] = 0;
		}


	}else if(strcasecmp(problem,"CONE") == 0){ 
	
		NDOF = 1;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]-10.0)<=1e-15 || fabs(coord[i][0])<=1e-15 || fabs(coord[i][1]-10.0)<=1e-15 ||fabs(coord[i][1])<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}
	}else if(strcasecmp(problem,"CONE2") == 0){ 
	
		NDOF = 1;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]+5.0)<=1e-15 || fabs(coord[i][0]-5.0)<=1e-15 || fabs(coord[i][1]+5.0)<=1e-15 ||fabs(coord[i][1]-5.0)<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}


	}else if(strcasecmp(problem,"CARTOLA") == 0){ 
	
		NDOF = 1;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]-1.0)<=1e-15 || fabs(coord[i][0])<=1e-15 || fabs(coord[i][1]-1.0)<=1e-15 ||fabs(coord[i][1])<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}

	} else if (strcasecmp(problem,"HEMKER")==0){
	
		NDOF = 1;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]*coord[i][0] + coord[i][1]*coord[i][1] - 1.0)<=1e-15 || fabs(coord[i][0] + 3.0)<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}

	}else if(strcasecmp(problem,"CAVITY") == 0){ 
		
		NDOF = 3;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]-1.0)<=1e-15 || fabs(coord[i][0])<=1e-15 || fabs(coord[i][1])<=1e-15 ||fabs(coord[i][1]-1.5)<=1e-15)
				mark[i] = 12;
			else if (fabs(coord[i][0]-0.5)<=1e-15 && fabs(coord[i][1]-0.5)<=1e-15)
				mark[i] = 2;
			else
				mark[i] = 0;
		}
	}else if(strcasecmp(problem,"CHANNEL") == 0){ 
		
		NDOF = 3;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0])<=1e-15 || fabs(coord[i][0]-30.0)<=1e-15 || fabs(coord[i][1])<=1e-15 ||fabs(coord[i][1]-1.5)<=1e-15
		           || (fabs(coord[i][0]-7.5)<=1e-15 && coord[i][1]<=0.75) || (fabs(coord[i][1]-0.75)<=1e-15 && coord[i][0]<=7.5))
				mark[i] = 12;
			else
				mark[i] = 0;
		}


	}else{
		printf("Problem not defined!\n");
		exit(1);
	}

	return NDOF;	
}


