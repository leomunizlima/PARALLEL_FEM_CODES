#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Build_K_F_DD(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int J, E, J1, J2, J3;
	int neq, nnodes, nel;
	double kx, ky, Be_x, Be_y, gamma, TwoA, invArea, Area, h_shock, Eu;
	double third=1.0/3.0, sixth = 1.0/6.0, seventh=1.0/7.0, twelfth= 1.0/12.0, three_twentieth=3./20.0, nine_twentieth=9./20.0, nine_forty_cents=9./40.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Beta[2], ke[3][3], XB, YB;
	double ue1, ue2, ue3, ueb, fe1, fe2, fe3, feb, Fe1, Fe2, Fe3;
	double KDD12, KDD13, KDD23, D12, D13, D23, C12, C13, C23, R11, R12, a1, a2, a3, invKBB, Khh[3][3], KhB[3], KBh[3], KBB, Fh[3], FB;
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

		XB = third*(X[0]+X[1]+X[2]);
		YB = third*(Y[0]+Y[1]+Y[2]);

		//Beta velocity calculation
		FemFunctions->Velocity(XB, YB, Beta);

		Be_x = Beta[0];
		Be_y = Beta[1];
		
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
	
		//Matrix Khh
		a1 = Be_x*y23 + Be_y*x32;
		a2 = Be_x*y31 + Be_y*x13;
		a3 = Be_x*y12 + Be_y*x21;

		C12 = y23*y31 + x32*x13;
        	C13 = y23*y12 + x32*x21;
		C23 = y31*y12 + x13*x21;

		D12 = 0.25*invArea*((kx * y23 * y31) + (ky * x32 * x13));
		D13 = 0.25*invArea*((kx * y23 * y12) + (ky * x32 * x21));
		D23 = 0.25*invArea*((kx * y31 * y12) + (ky * x13 * x21));

		R11 = sixth*gamma*Area;
		R12 = twelfth*gamma*Area;

		//----------------------  local length element --------------------------------------------------------/
		h_shock = FemFunctions->h_shock_calculation(Be_x, Be_y, ue1, ue2, ue3, y23, y31, y12, x32, x13, x21, Area);
		// ----------------------------------------------------------------------------------------------------/
		
		//---------------------- Shock capture parameter calculations ----------------------------------------------------------------/
		Eu = FemFunctions->ShockCapture(Re_old, delta_old, E, kx, ky, Be_x, Be_y, gamma, ue1, ue2, ue3, ueb, feb, y23, y31,y12, x32,x13,x21,invArea, h_shock); 
		//---------------------------------------------------------------------------------------------------------------------------/

		KDD12 = 0.25*Eu*invArea*C12;
		KDD13 = 0.25*Eu*invArea*C13;
		KDD23 = 0.25*Eu*invArea*C23;


		Khh[0][0] = -(D12 + D13) + sixth*a1 + R11 -(KDD12 + KDD13);
		Khh[0][1] =          D12 + sixth*a2 + R12 + KDD12;
		Khh[0][2] = 		 D13 + sixth*a3 + R12 + KDD13;
		Khh[1][0] =          D12 + sixth*a1 + R12 + KDD12;
		Khh[1][1] = -(D12 + D23) + sixth*a2 + R11 -(KDD12 + KDD23);
		Khh[1][2] = 		 D23 + sixth*a3 + R12 + KDD23;
		Khh[2][0] = 		 D13 + sixth*a1 + R12 + KDD13;
		Khh[2][1] = 		 D23 + sixth*a2 + R12 + KDD23;
		Khh[2][2] = -(D13 + D23) + sixth*a3 + R11 -(KDD13 + KDD23) ;

		// Matrix KhB
		KhB[0] = nine_twentieth*(0.5*a1 + a2 + a3) + three_twentieth*gamma*Area;
		KhB[1] = nine_twentieth*(a1 + 0.5*a2 + a3) + three_twentieth*gamma*Area;
		KhB[2] = nine_twentieth*(a1 + a2 + 0.5*a3) + three_twentieth*gamma*Area;


		// Matrix KBh
		KBh[0] =  nine_forty_cents*a1 + three_twentieth*gamma*Area;
		KBh[1] =  nine_forty_cents*a2 + three_twentieth*gamma*Area;
		KBh[2] =  nine_forty_cents*a3 + three_twentieth*gamma*Area;


		// Matrix KBB
		//Eu =  YZBeta_ShockCapture2(kx, ky, Be_x, Be_y, gamma, ue1, ue2, ue3, ueb, feb, y23, y31,y12, x32,x13,x21,invArea, h_shock); 
		
		KBB = -9*nine_forty_cents*(D12+D13+D23) + seventh*three_twentieth*(a1+ a2+ a3) + 9*seventh*nine_forty_cents*Area - Eu*invArea*(C12 + C13 + C23)/360.;
		invKBB = 1.0/KBB;

		// Vector Fh
		Fh[0] = twelfth*Area*(2*fe1 + fe2 + fe3);
		Fh[1] = twelfth*Area*(fe1 + 2*fe2 + fe3);
		Fh[2] = twelfth*Area*(fe1 + fe2 + 2*fe3);

		// Vector FB
		FB = three_twentieth*Area*(fe1 + fe2 + fe3);

		ke[0][0] = Khh[0][0] - invKBB*KhB[0]*KBh[0];
		ke[0][1] = Khh[0][1] - invKBB*KhB[0]*KBh[1];
		ke[0][2] = Khh[0][2] - invKBB*KhB[0]*KBh[2];
		ke[1][0] = Khh[1][0] - invKBB*KhB[1]*KBh[0];
		ke[1][1] = Khh[1][1] - invKBB*KhB[1]*KBh[1];
		ke[1][2] = Khh[1][2] - invKBB*KhB[1]*KBh[2];
		ke[2][0] = Khh[2][0] - invKBB*KhB[2]*KBh[0];
		ke[2][1] = Khh[2][1] - invKBB*KhB[2]*KBh[1];
		ke[2][2] = Khh[2][2] - invKBB*KhB[2]*KBh[2];

		Fe1 = Fh[0] - invKBB*KhB[0]*FB;
		Fe2 = Fh[1] - invKBB*KhB[1]*FB;
		Fe3 = Fh[2] - invKBB*KhB[2]*FB;

		// Boundary conditions on vector F and Fill Matrix Diagonal 
		F_assembly(Node, J1, J2, J3, Fe1, Fe2, Fe3, F, ke, X, Y, FemFunctions);

		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, E, ke);

	}

	free(U);

	return 0;
}


