#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int setmark_according_boundary(int nnodes, double **coord, int *mark, char *problem);

int main(int argc, char **argv)
{
	char inclination[30], label[300], FileName[300], problem[100];
	int M, N;          // M: number of subdivisions in x e N: number of subdivisions in y
	FILE *File;
	int nosx, nosy, nnodes, nel, nodeCont;
	double xInf, yInf, xSup, ySup, hx, hy;
	int k, i, j, e, Xcria, Ycria;
	int no1, no2, no3;
	int print_bound[16][4]= {{1,1,1,1}/*mark 0 */, {1,1,1,0}/*mark 1 */, {1,1,0,1}/*mark 2 */, {1,1,0,0}/*mark 3 */,
				{1,0,1,1}/*mark 4 */, {1,0,1,0}/*mark 5 */, {1,0,0,1}/*mark 6 */, {1,0,0,0}/*mark 7 */,
				{0,1,1,1}/*mark 8 */, {0,1,1,0}/*mark 9 */, {0,1,0,1}/*mark 10*/, {0,1,0,0}/*mark 11*/,
				{0,0,1,1}/*mark 12*/, {0,0,1,0}/*mark 13*/, {0,0,0,1}/*mark 14*/, {0,0,0,0}/*mark 15*/};
	
	if (argc != 3){
		printf("Use ./meshGenerator <PROBLEM>.txt <PROBLEM>\n");
		exit(1);
	}
	
	strcpy(problem,argv[2]);
	
	// Input file reading
	File = fopen(argv[1], "r");
	if (File == NULL){
		printf("File %s not found!\n", argv[1]);
		exit(1);
	}
	fscanf(File, "%s\t:%[^\n]", inclination, label);
	fscanf(File, "%d\t:%[^\n]", &M, label);
	fscanf(File, "%d\t:%[^\n]", &N, label);
	fscanf(File, "%d\t:%[^\n]", &Xcria, label);
	fscanf(File, "%d\t:%[^\n]", &Ycria, label);
	fscanf(File, "%d\t:%[^\n]", &nodeCont, label);
	
	double *x, *y; // vectors coordinates of the nodes forming the boundary
	x = calloc(nodeCont, sizeof(double));
	y = calloc(nodeCont, sizeof(double));

	for(i = 0; i < nodeCont; i++){
		fscanf(File, "%lf\t%lf\t:%[^\n]", &x[i], &y[i], label);
	}
	fclose(File);
	

	printf("Inclination: %s\n", inclination);
	printf("M: %d\n", M);
	printf("N: %d\n", N);
	printf("Xcria: %d\n", Xcria);
	printf("Ycria: %d\n", Ycria);
	printf("NodeCont: %d\n", nodeCont);
	for(i = 0; i < nodeCont; i++){
		printf("x%d: %lf \t y%d: %lf\n", i, x[i], i, y[i]);
	}

	// Starting calculations
	nosx = M+1;
	nosy = N+1;
	nnodes = nosx*nosy;

	xInf = x[0];
	yInf = y[0];
	xSup = x[0];
	ySup = y[0];

	for(i = 0; i < nodeCont; i++){
		if(xInf > x[i]){
			xInf = x[i];
		}
		if(yInf > y[i]){
			yInf = y[i];
		}
	}	

	for(i = 0; i < nodeCont; i++){
		if(xSup < x[i]){
			xSup = x[i];
		}
		if(ySup < y[i]){
			ySup = y[i];
		}
	}	

	free(x);
	free(y);
	

	nel = M*N*2;
	hx = (xSup - xInf)/M;
	hy = (ySup - yInf)/N;
	
	double **coord;
	int **ien, *mark;

	coord = calloc(nnodes,sizeof(double*));
	ien = calloc(nel,sizeof(int*));
	mark = calloc(nnodes,sizeof(int));

	for (i=0; i<nnodes; i++)
		coord[i] = calloc(2,sizeof(double));
	
	for (i=0; i<nel; i++)
		ien[i] = calloc(3,sizeof(int));
	

	//Right-tilted mesh
	if (strcmp(inclination,"right") == 0){
		no1 = 0;
		no2 = 0;
		no3 = 0;

		j = 0;
		k = (2 * nosx) - 2;         
		for (e = 0; e < nel; e++){
			j++;       
			if (e % 2 == 0){       
				no1++;
				no2 = no1 + 1;
				no3 = no2 + nosx;
			}else{                  
				no2 = no3;
				no3 = no2 - 1;
			}
			ien[e][0] = no1-1;
			ien[e][1] = no2-1;
			ien[e][2] = no3-1;

			if (j == k){
				no1++;
				j = 0;
			}                                
		}
	//Left-tilted mesh
	}else{
		no1 = 1;
		no2 = 0;
		no3 = 0;

		j = 0;
		k = (2 * nosx) - 2;         
		for (e = 0; e < nel; e++){
			j++;  
			if (e % 2 == 0){       
				no1 = no1;
				no2 = no1 + 1;
				no3 = no1 + nosx;
			}else{      
				no1 = no2;
				no2 = no3 + 1;
				no3 = no2 - 1;
			}   
			ien[e][0] = no1-1;
			ien[e][1] = no2-1;
			ien[e][2] = no3-1;   

			if (j == k){
				no1++;
				j = 0;
			}                         
		}
	}

	// coordenators in global nodes
	k = 0;
	for (j = 0; j < nosy; j++){
		for (i = 0; i < nosx; i++){
			coord[k][0] = xInf + i*hx;
			coord[k][1] = yInf + j*hy;
			k++;
		}
	}

	//setting mark according to boundary conditions
	int NDOF;
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
			fprintf(File, "%d \t %d \t %d\n", ien[i][0], ien[i][1], ien[i][2]);
	}
	else {
		for(i = 0; i < nel; i++)
			fprintf(File, "%d \t %d \t %d \t -1\n", ien[i][0], ien[i][1], ien[i][2]);
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
	
}//end main


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
	}else if(strcasecmp(problem,"SOD") == 0){
	
		double max;
	
		NDOF = 4;
		max = coord[0][1];
		for (i=0; i<nnodes;i++){
			if (coord[i][1]>max)
				max = coord[i][1];	
		}
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0])<=1e-15 || fabs(coord[i][0]-1.0)<=1e-15)
				mark[i] = 15;
			else if (fabs(coord[i][1])<=1e-15 || fabs(coord[i][1]-max)<=1e-15)
				mark[i] = 2;
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

	}else if(strcasecmp(problem,"PUDIM") == 0){ 

		NDOF = 1;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]-1.0)<=1e-15 || fabs(coord[i][0]+1.0)<=1e-15 || fabs(coord[i][1]-1.0)<=1e-15 ||fabs(coord[i][1]+1.0)<=1e-15 ||
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

	}else if(strcasecmp(problem,"CARTOLA") == 0){ 
	
		NDOF = 1;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]-1.0)<=1e-15 || fabs(coord[i][0])<=1e-15 || fabs(coord[i][1]-1.0)<=1e-15 ||fabs(coord[i][1])<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}
	}else if(strcasecmp(problem,"TESTE") == 0){ 
	
		NDOF = 1;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]-1.0)<=1e-15 || fabs(coord[i][0])<=1e-15 || fabs(coord[i][1]-1.0)<=1e-15 ||fabs(coord[i][1])<=1e-15)
				mark[i] = 15;
			else
				mark[i] = 0;
		}
	}else if(strcasecmp(problem,"CAVITY") == 0){ 
		
		NDOF = 3;
		for(i=0;i<nnodes;i++){
			if (fabs(coord[i][0]-1.0)<=1e-15 || fabs(coord[i][0])<=1e-15 || fabs(coord[i][1]-1.0)<=1e-15 ||fabs(coord[i][1])<=1e-15)
				mark[i] = 12;
			else if (fabs(coord[i][0]-0.5)<=1e-15 && fabs(coord[i][1]-0.5)<=1e-15)
				mark[i] = 2;
			else
				mark[i] = 0;
		}

	}else{
		printf("Problem not defined!\n");
		exit(1);
	}

	return NDOF;	
}





