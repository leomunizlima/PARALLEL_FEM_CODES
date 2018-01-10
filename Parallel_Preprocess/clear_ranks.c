#include<stdlib.h>
#include<stdio.h>

int main(int argc, char **argv)
{
	if (argc<2){
		printf("Use ./clear_ranks <path of directories> <list of directories to be cleaned>\n");
		exit(1);
	}

	int i, number_of_directories = argc-2;
	char command[300];
	
	for (i=0;i<number_of_directories;i++){
		command[0]='\0';
		sprintf(command,"rm -r %s/%02d/RANK_*.dat\n",argv[1],atoi(argv[i+2]));
		system(command);
	}
	
	return 0;

}
