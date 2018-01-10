# include "explosion.h"

int EXPLOSION_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u){

	int I, I1, I2, I3, I4, nnodes;
	double x, y;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++){

		x = Node[I].x;	
		y = Node[I].y;	

		I1 = Node[I].invP_id[0];
		I2 = Node[I].invP_id[1];
		I3 = Node[I].invP_id[2];
		I4 = Node[I].invP_id[3];

		//Sphere in 
		if ((x-1)*(x-1) + (y-1)*(y-1)< 0.16){//<0.14
			if (Node[I].Type[0] != 0)
				u[I1] = 1.0; 
			if (Node[I].Type[1] != 0)
				u[I2] = 0.0; 
			if (Node[I].Type[2] != 0)
				u[I3] = 0.0; 
			if (Node[I].Type[3] != 0)
				u[I4] = 2.5; // rho*e 
		}
		//Sphere out
		else{
			if (Node[I].Type[0] != 0)
				u[I1] = 0.125; 
			if (Node[I].Type[1] != 0)
				u[I2] = 0.0; 
			if (Node[I].Type[2] != 0)
				u[I3] = 0.0; 
			if (Node[I].Type[3] != 0)
				u[I4] = 0.25;  //rho * e

		}	
	}

	return 0;
}


