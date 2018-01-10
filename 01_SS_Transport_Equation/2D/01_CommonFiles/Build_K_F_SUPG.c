#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Build_K_F_SUPG(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int J, E, J1, J2, J3;
	int neq, nnodes, nel;
	double kx, ky, Be_x, Be_y, gamma, TwoA, invArea, Area, norm, Be2_x, Be2_y, tau, h, h_shock, Eu;
	double third=1.0/3.0, sixth = 1.0/6.0, twelfth= 1.0/12.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Beta[2], ke[3][3];
	double ue1, ue2, ue3, ueb, fe1, fe2, fe3, feb, Fe1, Fe2, Fe3, Rpg11, Rpg21, Rpg31, Cpg12, Cpg13, Cpg23, Csc12, Csc13, Csc23;
	double D12, D13, D23, C11, C12, C13, R11, R12, XB, YB;
	double *F, *u, *U; // U represent values of unknows and boundary conditions. 
	double *delta_old = FemStructs->delta_old;
	double *Re_old = FemStructs->Re_old;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	
	u = FemStructs->u;
	F = FemStructs->f;
	nel = Parameters->nel;
	neq = Parameters->neq;
	nnodes = Parameters->nnodes;
	kx = FemFunctions->Condutivity();
	ky = FemFunctions->Condutivity();
	gamma = FemFunctions->Reaction();
	setzeros(Parameters,MatrixData);
//=============================================================================================================================
/*	FILE *Out;
	char Filename[200];
	sprintf(Filename,"Matrices_%d_%s.dat",Parameters->RANK, Parameters->reordering);
	Out=fopen(Filename, "w");*/
//===================================================================================================================================

	for (J=0;J<neq;J++)
		F[J] = 0;

	U = mycalloc("U of Build_K_F_SUPG",nnodes,sizeof(double));	

	for (J=0; J<nnodes; J++){
		if (Node[J].Type==0)
			U[J] = FemFunctions->upresc(Node[J].x,Node[J].y);
		else
			U[J] = u[Node[J].invP_id];
	}

	for (E=0; E<nel; E++)
	{
		J1 = Element[E].Vertex[0];
		J2 = Element[E].Vertex[1];
		J3 = Element[E].Vertex[2];

		X[0] = Node[J1].x;
		X[1] = Node[J2].x;
		X[2] = Node[J3].x;
		Y[0] = Node[J1].y;
		Y[1] = Node[J2].y;
		Y[2] = Node[J3].y;

		y23 = Y[1] - Y[2];
		y31 = Y[2] - Y[0];
		y12 = Y[0] - Y[1];

		x32 = X[2] - X[1];
		x13 = X[0] - X[2];
		x21 = X[1] - X[0];

		TwoA = x21*y31 - x13*y12;
		Area = 0.5*TwoA;
		invArea = 1.0/Area;
		h = sqrt(fabs(TwoA));

		XB = third*(X[0]+X[1]+X[2]);
		YB = third*(Y[0]+Y[1]+Y[2]);

		//Beta velocity calculation
		FemFunctions->Velocity(XB, YB, Beta);

		Be_x = Beta[0];
		Be_y = Beta[1];
		Be2_x = Be_x*Be_x;
		Be2_y = Be_y*Be_y;

		norm = sqrt(Be2_x + Be2_y);

		if (norm>0)
			tau = 0.5*h/norm;
		else
			tau = 0;

		//Calculation of Font/////////
		fe1 = FemFunctions->Font(X[0], Y[0], kx, gamma, Be_x, Be_y);
		fe2 = FemFunctions->Font(X[1], Y[1], kx, gamma, Be_x, Be_y);
		fe3 = FemFunctions->Font(X[2], Y[2], kx, gamma, Be_x, Be_y);

		feb = third*(fe1+fe2+fe3);

		//Calculation of ue in element e 
		ue1 = U[J1];
		ue2 = U[J2];
		ue3 = U[J3];

		ueb = third * (ue1 + ue2 + ue3);

		// Galerkin diffusion matrix 
		D12 = 0.25*invArea*((kx * y23 * y31) + (ky * x32 * x13));
		D13 = 0.25*invArea*((kx * y23 * y12) + (ky * x32 * x21));
		D23 = 0.25*invArea*((kx * y31 * y12) + (ky * x13 * x21));

		// Galerkin Convection matrix
		C11 = sixth*(Be_x*y23 + Be_y*x32);
		C12 = sixth*(Be_x*y31 + Be_y*x13);
		C13 = sixth*(Be_x*y12 + Be_y*x21);

		// Galerkin reaction matrix
		R11 = sixth*(gamma*Area);
		R12 = twelfth*(gamma * Area);

		// SUPG convection matrix
	  	Cpg12 = 0.25*invArea*tau * (y23*(Be_x*Be_x*y31 + Be_x*Be_y*x13) + x32*(Be_x*Be_y*y31 + Be_y*Be_y*x13));
	   	Cpg13 = 0.25*invArea*tau * (y23*(Be_x*Be_x*y12 + Be_x*Be_y*x21) + x32*(Be_x*Be_y*y12 + Be_y*Be_y*x21));
	   	Cpg23 = 0.25*invArea*tau * (y31*(Be_x*Be_x*y12 + Be_x*Be_y*x21) + x13*(Be_x*Be_y*y12 + Be_y*Be_y*x21));

	       // SUPG reaction matrix
	   	Rpg11 = sixth*tau*gamma * (Be_x*y23 + Be_y*x32);
	   	Rpg21 = sixth*tau*gamma * (Be_x*y31 + Be_y*x13);
	   	Rpg31 = sixth*tau*gamma * (Be_x*y12 + Be_y*x21);

	   	
		//----------------------  local length element --------------------------------------------------------/
		h_shock = FemFunctions->h_shock_calculation(Be_x, Be_y, ue1, ue2, ue3, y23, y31, y12, x32, x13, x21, Area);
		// ----------------------------------------------------------------------------------------------------/

		//---------------------- Shock capture parameter calculations ----------------------------------------------------------------/
		Eu = FemFunctions->ShockCapture(Re_old, delta_old, E, kx, ky, Be_x, Be_y, gamma, ue1, ue2, ue3, ueb, feb, y23, y31,y12, x32,x13,x21,invArea, h_shock); 
		//---------------------------------------------------------------------------------------------------------------------------/

		// Shock capture correction matrix
		Csc12 = Eu*0.25*invArea*(y23*y31 + x32*x13);
	   	Csc13 = Eu*0.25*invArea*(y23*y12 + x32*x21);
	   	Csc23 = Eu*0.25*invArea*(y31*y12 + x13*x21);

	    	// Matrix ke
		ke[0][0] = -(D12 + D13) + C11 - (Cpg12 + Cpg13) + R11 + Rpg11 - (Csc12 + Csc13);
		ke[0][1] = 	    D12 + C12 + 	Cpg12 + R12 + Rpg11 + Csc12;
		ke[0][2] = 		D13 + C13 + 	Cpg13 + R12 + Rpg11 + Csc13;
		ke[1][0] = 		D12 + C11 + 	Cpg12 + R12 + Rpg21 + Csc12;
		ke[1][1] = -(D12 + D23) + C12 - (Cpg12 + Cpg23) + R11 + Rpg21 - (Csc12 + Csc23);
		ke[1][2] = 		D23 + C13 + Cpg23 + R12 + Rpg21 + Csc23;
		ke[2][0] = 		D13 + C11 + Cpg13 + R12 + Rpg31 + Csc13;
		ke[2][1] = 		D23 + C12 + Cpg23 + R12 + Rpg31 + Csc23;
		ke[2][2] = -(D13 + D23) + C13 - (Cpg13 + Cpg23) + R11 + Rpg31 - (Csc13 + Csc23);

/*		int i,j;
//fprintf(Out,"E=%d\n x1=%.20lf x2=%.20lf x3=%.20lf y1=%.20lf y2=%.20lf y3=%.20lf\n",E,X[0],X[1],X[2],Y[0],Y[1],Y[2]); 
		fprintf(Out,"E=%d\n",E); 
		for (i=0; i<3; i++)
			for (j=0; j<3; j++)
				fprintf(Out,"ke[%d][%d]=%.20lf\n",i,j,ke[i][j]);
*/				
			
		

		
		// F contribution
		Fe1 = Area*twelfth*(2*fe1 +  fe2  +  fe3 ) + sixth*tau*(fe1+fe2+fe3)*(Be_x*y23 + Be_y*x32);
		Fe2 = Area*twelfth*(  fe1 + 2*fe2 +  fe3 ) + sixth*tau*(fe1+fe2+fe3)*(Be_x*y31 + Be_y*x13);
		Fe3 = Area*twelfth*(  fe1 +  fe2  + 2*fe3) + sixth*tau*(fe1+fe2+fe3)*(Be_x*y12 + Be_y*x21);


		// Boundary conditions on vector F and Fill Matrix Diagonal 
		F_assembly(Node, J1, J2, J3, Fe1, Fe2, Fe3, F, ke, X, Y, FemFunctions);
		
		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, E, ke);
		
	}

/*	int I;
	for (I=0; I<Parameters->nnodes; I++)
		if (Node[I].Type==1)
			fprintf(Out,"F[%d]=%.20lf\n",Node[I].invP_id,F[Node[I].invP_id]);
	fclose(Out);	
*/
	free(U);
	return 0;
}



