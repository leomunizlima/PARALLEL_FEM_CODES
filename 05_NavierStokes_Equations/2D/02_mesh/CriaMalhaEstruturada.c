#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//========= ALOCACAO DINAMICA DE MATRIZ DOUBLE ========
double **alocMatDouble(int Linhas,int Colunas){ 					//Recebe a quantidade de Linhas e Colunas como Parâmetro
	int i,j; 													//Variáveis Auxiliares
	double **m = (double**)malloc(Linhas*sizeof(double*)); 	//Aloca um Vetor de Ponteiros
	for (i = 0; i < Linhas; i++){ 								//Percorre as linhas do Vetor de Ponteiros
		m[i] = (double*)malloc(Colunas*sizeof(double));			//Aloca um Vetor de Inteiros para cada posição do Vetor de Ponteiros.
		for (j = 0; j < Colunas; j++){ 							//Percorre o Vetor de Inteiros atual.
			//m[i][j] = ((j%2 == 0 )? 0.: 1.);		//LINHA TESTE, USO FA FUNCAO MODULO
			m[i][j] = 0.; 										//Inicializa com 0.
		}
	}
return m;                                          				 //Retorna o Ponteiro para a Matriz Alocada
}
//======================================================
//========= ALOCACAO DINAMICA DE MATRIZ INT ========
int **alocMatInt(int Linhas,int Colunas){ 					//Recebe a quantidade de Linhas e Colunas como Parâmetro
	int i,j; 													//Variáveis Auxiliares
	int **m = (int**)malloc(Linhas*sizeof(int*)); 	//Aloca um Vetor de Ponteiros
	for (i = 0; i < Linhas; i++){ 								//Percorre as linhas do Vetor de Ponteiros
		m[i] = (int*)malloc(Colunas*sizeof(int));			//Aloca um Vetor de Inteiros para cada posição do Vetor de Ponteiros.
		for (j = 0; j < Colunas; j++){ 							//Percorre o Vetor de Inteiros atual.
			//m[i][j] = ((j%2 == 0 )? 0.: 1.);		//LINHA TESTE, USO FA FUNCAO MODULO
			m[i][j] = 0.; 										//Inicializa com 0.
		}
	}
return m;                                          				 //Retorna o Ponteiro para a Matriz Alocada
}
//===================================================
//========= ALOCACAO DINAMICA DE VETOR INT===============
int *alocVetInt(int Comp){ 					//Recebe a quantidade de Linhas e Colunas como Parâmetro
	int i; 													//Variáveis Auxiliares
	int *vet;
	vet = (int*)malloc(Comp*sizeof(int)); 	//Aloca um Vetor de Ponteiros
	for (i = 0; i < Comp; i++) 								//Percorre as linhas do Vetor de Ponteiros
		vet[i] = 0; 										//Inicializa com 0.
return vet;                          		 //Retorna o Ponteiro para a Matriz Alocada
}

//============PRINCIPAL===========================
int main(int argc, char **argv)
{
	char inclination[30], label[300], FileName[300], problem[20], name[25];
	int M, N;          // M: number of subdivisions in x e N: number of subdivisions in y
	int markNode[6];   // vector that holds the marks of the nodes , the node count begins in the lower left and runs counter-clockwise
	int markSide[6];   // vector that holds the marks of the contours , the contour of the count begins at the bottom and rotates counter- clockwise
	double x[6], y[6]; // vectors coordinates of the nodes forming the boundary
	FILE *File;
	int nosx, nosy, nnodes, nel, nodeCont, edgeCont;
	double xInf, yInf, xSup, ySup, hx, hy, Y, X, cy, cx;
	int k, i, j, e, Xcria, Ycria;
	int no1, no2, no3, size;
	
	if (argc != 3){
		printf("Use ./CriaEstruturada EntradaMalhaEstruturada.txt <PROBLEM>\n");
		exit(1);
	}
	
	strcpy(problem,argv[2]);
	
	// Leitura de arquivo de entrada
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
	for(i = 0; i < nodeCont; i++){
		fscanf(File, "%lf\t%lf\t%d\t:%[^\n]", &x[i], &y[i], &markNode[i], label);
	}
	fscanf(File, "%d\t:%[^\n]", &edgeCont, label);
	for(i = 0; i < edgeCont; i++){
		fscanf(File, "%d\t:%[^\n]", &markSide[i], label);
	}
	fclose(File);
	
	// Verificando se leu certo
	printf("Inclination: %s\n", inclination);
	printf("M: %d\n", M);
	printf("N: %d\n", N);
	printf("Xcria: %d\n", Xcria);
	printf("Ycria: %d\n", Ycria);
	printf("NodeCont: %d\n", nodeCont);
	for(i = 0; i < nodeCont; i++){
		printf("x%d: %lf \t y%d: %lf \t Mark: %d\n", i, x[i], i, y[i], markNode[i]);
	}
	printf("EdgeCont: %d\n", edgeCont);
	for(i = 0; i < edgeCont; i++){
		printf("Mark: %d\n", markSide[i]);
	}

	// Iniciando os calculos
	nosx = M+1;
	nosy = N+1;
	nnodes = nosx*nosy;

	xInf = x[0];
	yInf = y[0];
	xSup = x[0];
	ySup = y[0];
	for(i = 0; i < nodeCont; i++){
		if(xSup < x[i]){
			xSup = x[i];
		}
		if(ySup < y[i]){
			ySup = y[i];
		}
	}	

	nel = M*N*2;
	hx = (xSup - xInf)/M;
	hy = (ySup - yInf)/N;
	
	//double coord[nnodes][2];
	double **coord;
	coord = alocMatDouble(nnodes,2);
	//int ien[nel][3], mark[nnodes];
	int **ien, *mark;
	ien = alocMatInt(nel,3);
	mark = alocVetInt(nnodes);

	/*for(i = 0; i < nnodes; i++){
		mark[i] = 0;
	}*/

	//malha inclinada para direita
	if (strcmp(inclination,"direita") == 0){
		no1 = 0;
		no2 = 0;
		no3 = 0;

		j = 0;
		k = (2 * nosx) - 2;         // numero de elementos em x
		for (e = 1; e <= nel; e++){
			j++;       
			if (e % 2 == 1){        // elemento de numero impar
				no1++;
				no2 = no1 + 1;
				no3 = no2 + nosx;
			}else{                  // elemento de numero par
				no2 = no3;
				no3 = no2 - 1;
			}
			ien[e - 1][0] = no1-1;
			ien[e - 1][1] = no2-1;
			ien[e - 1][2] = no3-1;

			if (j == k){
				no1++;
				j = 0;
			}                                
		}
	}else{  // malha inclinada para esquerda
		no1 = 1;
		no2 = 0;
		no3 = 0;

		j = 0;
		k = (2 * nosx) - 2;         // numero de elementos em x
		for (e = 1; e <= nel; e++){
			j++;  
			if (e % 2 == 1){       // elemento de numero impar
				no1 = no1;
				no2 = no1 + 1;
				no3 = no1 + nosx - 1;
			}else{      // elemento de numero par
				no1 = no2;
				no2 = no3 + 1;
				no3 = no2 - 1;
			}   
			ien[e-1][0] = no1-1;
			ien[e-1][1] = no2-1;
			ien[e-1][2] = no3-1;   

			if (j == k){
				no1++;
				j = 0;
			}                         
		}
	}

	// coordenadas dos nos globais
	k = 0;
	for (j = 0; j < nosy; j++){
		for (i = 0; i < nosx; i++){
			coord[k][0] = xInf + i*hx;
			coord[k][1] = yInf + j*hy;
			k++;
		}
	}

	k = 0;
	Y = yInf; X = xInf;
	mark[k] = markNode[0];
	k++;
	
	for (X += hx; X < x[1]; X += hx){
		mark[k] = markSide[0];
		k++;
	}
	if (Xcria==3) {
		mark[k] = markNode[1];
		X += hx;
		k++;
	}
	for (; X < xSup; X += hx, k++){
		mark[k] = markSide[1];
	}
	
	mark[k] = markNode[Xcria-1];
	k++;
	
	for (Y += hy; Y < ySup; Y += hy, k++){
		mark[k] = markSide[Xcria*Ycria-1];
		k+=M;
		mark[k] = markSide[Xcria-1];
	}
	
	Y = ySup; X = xInf;
	mark[k] = markNode[Xcria*Ycria-1];
	k++;
	
	for (X += hx; X < x[1]; X += hx){
		mark[k] = markSide[Xcria*Ycria-2];
		k++;
	}
	
	if (Xcria==3) {
		mark[k] = markNode[Xcria*Ycria-2];
		X += hx;
		k++;
	}
	for (; X < xSup; X += hx){
		mark[k] = markSide[3];
		k++;
	}
	
	mark[k] = markNode[Xcria];

	// Gerando arquivo de saida
	size = strlen(argv[1]);
	strncpy(name, argv[1], size-10);
	name[size-10] = '\0';
	sprintf(FileName,"%s_%d_%d.dat", name, nnodes, nel);
	File = fopen(FileName, "w");
	if (File == NULL)
	{
		printf("File %s not found!\n", FileName);
		exit(1);
	}

	if(strcasecmp(problem,"OBLIQUO") == 0){
		fprintf(File, "%d\n", nnodes);
		for(i = 0; i < nnodes; i++){
			if(mark[i] == 1){
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d \t %d\n", coord[i][0], coord[i][1], 1, 1, 0, 1);
			}else if(mark[i] == 3){
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d \t %d\n", coord[i][0], coord[i][1], 0, 0, 0, 0);				
			}else{
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d \t %d\n", coord[i][0], coord[i][1], 1, 1, 1, 1);
			}
		}
		fprintf(File, "%d\n", nel);
		for(i = 0; i < nel; i++){
			fprintf(File, "%d \t %d \t %d \t -1\n", ien[i][0], ien[i][1], ien[i][2]);
		}
	}else if(strcasecmp(problem,"REFLETIDO") == 0){
		fprintf(File, "%d\n", nnodes);
		for(i = 0; i < nnodes; i++){
			if(mark[i] == 1){
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d \t %d\n", coord[i][0], coord[i][1], 1, 1, 0, 1);
			}else if(mark[i] == 3){
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d \t %d\n", coord[i][0], coord[i][1], 0, 0, 0, 0);				
			}else{
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d \t %d\n", coord[i][0], coord[i][1], 1, 1, 1, 1);
			}
		}
		fprintf(File, "%d\n", nel);
		for(i = 0; i < nel; i++){
			fprintf(File, "%d \t %d \t %d \t -1\n", ien[i][0], ien[i][1], ien[i][2]);
		}
	}else if(strcasecmp(problem,"SOD") == 0){
		fprintf(File, "%d\n", nnodes);
		for(i = 0; i < nnodes; i++){
			if(mark[i] == 1){
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d \t %d\n", coord[i][0], coord[i][1], 1, 1, 0, 1);
			}else if(mark[i] == 2){
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d \t %d\n", coord[i][0], coord[i][1], 0, 0, 0, 0);				
			}else{
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d \t %d\n", coord[i][0], coord[i][1], 1, 1, 1, 1);
			}
		}
		fprintf(File, "%d\n", nel);
		for(i = 0; i < nel; i++)
			fprintf(File, "%d \t %d \t %d \t -1\n", ien[i][0], ien[i][1], ien[i][2]);
	}else if(strcasecmp(problem,"CAVITY") == 0){
		fprintf(File, "%d\n", nnodes);
		for(i = 0; i < nnodes; i++){
			cx = coord[i][0];
			cy = coord[i][1];
			if(mark[i] == 1 || mark[i] == 2 || mark[i] == 3 || mark[i] == 4){
				if((fabs(cx-0.5) <= 1e-15) && (fabs(cy-0.0) <= 1e-15))
					fprintf(File, "%lf \t %lf \t %d \t %d \t %d\n", cx, cy, 0, 0, 0);
				else
					fprintf(File, "%lf \t %lf \t %d \t %d \t %d\n", cx, cy, 0, 0, 1);
			}else{
					fprintf(File, "%lf \t %lf \t %d \t %d \t %d\n", cx, cy, 1, 1, 1);
			}
		}
		fprintf(File, "%d\n", nel);
		for(i = 0; i < nel; i++){
			fprintf(File, "%d \t %d \t %d \t -1\n", ien[i][0], ien[i][1], ien[i][2]);
		}
	}else if(strcasecmp(problem,"EXATA") == 0){
		fprintf(File, "%d\n", nnodes);
		for(i = 0; i < nnodes; i++){
			cx = coord[i][0];
			cy = coord[i][1];
			if(mark[i] == 1 || mark[i] == 2 || mark[i] == 3 || mark[i] == 4)			
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d\n", cx, cy, 0, 0, 0);
			else
				fprintf(File, "%lf \t %lf \t %d \t %d \t %d\n", cx, cy, 1, 1, 1);
		}
		fprintf(File, "%d\n", nel);
		for(i = 0; i < nel; i++){
			fprintf(File, "%d \t %d \t %d \t -1\n", ien[i][0], ien[i][1], ien[i][2]);
		}
	}

	fclose(File);
	
	return 0;
	
}//end main
