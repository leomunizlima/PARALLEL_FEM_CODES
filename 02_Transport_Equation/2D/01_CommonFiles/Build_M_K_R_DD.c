#include "TranspEquation.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_M_K_R_DD(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int i, j, E, J1, J2, J3;
	int nel, NEQ;
	double kx, ky, Be_x, Be_y, gamma, TwoA, invArea, Area, h_shock, Eu, alpha, dt;
	double w=0.5;
	double third=1.0/3.0, sixth = 1.0/6.0, seventh=1.0/7.0, twelfth= 1.0/12.0,three_twentieth=3./20.0, nine_twentieth=9./20.0, nine_forty_cents=9./40.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Beta[2], Ae[3][3], ke[3][3], me[3][3], XB, YB;
	double ue[3], ueb, due[3], dueb, fe1, fe2, fe3, feb, *U, *dU;
	double KDD12, KDD13, KDD23, D12, D13, D23, C12, C13, C23, R11, R12, a1, a2, a3, invKBB, Mhh_1_12, MhB_3_20, MBh_3_20,
	Khh[3][3], KhB[3], KBh[3], KBB, Fh[3], FB, cbtil, cb;
	double *R = FemStructs->F;
	double *cbtilo = FemStructs->delta_old;
	int **lm = FemStructs->lm;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	nel = Parameters->nel;
	NEQ = Parameters->NEQ;
	kx = FemFunctions->Condutivity();
	ky = FemFunctions->Condutivity();
	gamma = FemFunctions->Reaction();

	alpha = Parameters->Alpha_Build;
	dt = Parameters->DeltaT_Build;

	dzero(NEQ+1, R);
	setzeros(Parameters,MatrixData);

	U = (double*) mycalloc("U of 'Build_M_K_R_DD'", Parameters->nnodes, sizeof(double));
	dU = (double*) mycalloc("dU of 'Build_M_K_R_DD'", Parameters->nnodes, sizeof(double));
	eval_U_dU(Parameters,FemStructs,FemFunctions,U,dU);

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

		//Calculation of u and du in element e
		ue[0] = U[J1];
		ue[1] = U[J2];
		ue[2] = U[J3];
		due[0] = dU[J1];
		due[1] = dU[J2];
		due[2] = dU[J3];
		
		ueb = third * (ue[0] + ue[1] + ue[2]);
	   	dueb = third * (due[0] + due[1] + due[2]);

		//Matrix Mhh
		Mhh_1_12 = twelfth*Area;

		//Matrix MhB
		MhB_3_20 = three_twentieth*Area;

		//Matrix MBh
		MBh_3_20 = MhB_3_20;


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
	
		/***************** local length element **********************/
		h_shock = FemFunctions->h_shock(Be_x, Be_y, ue[0], ue[1], ue[2], y23, y31, y12, x32, x13, x21, Area);
		/**************************************************************/

		/**************************** Shock capture parameter calculations***************************************/
		cbtil = FemFunctions->ShockCapture(kx, ky, Be_x, Be_y, gamma, ue[0], ue[1], ue[2], ueb, dueb, feb, y23, y31,y12, x32,x13,x21,invArea, h_shock); 
		cb = w*cbtil + (1-w)*cbtilo[E];
		cbtilo[E] = cbtil;
		Eu = cb;
		/********************************************************************************************************/

		KDD12 = 0.25*Eu*invArea*C12;
		KDD13 = 0.25*Eu*invArea*C13;
		KDD23 = 0.25*Eu*invArea*C23;

		Khh[0][0] = -(D12 + D13) + sixth*a1 + R11 -(KDD12 + KDD13);
		Khh[0][1] =          D12 + sixth*a2 + R12 + KDD12;
		Khh[0][2] = 		D13 + sixth*a3 + R12 + KDD13;
		Khh[1][0] =          D12 + sixth*a1 + R12 + KDD12;
		Khh[1][1] = -(D12 + D23) + sixth*a2 + R11 -(KDD12 + KDD23);
		Khh[1][2] = 	     D23 + sixth*a3 + R12 + KDD23;
		Khh[2][0] = 	     D13 + sixth*a1 + R12 + KDD13;
		Khh[2][1] = 	     D23 + sixth*a2 + R12 + KDD23;
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
		KBB = -9*nine_forty_cents*(D12+D13+D23) + seventh*three_twentieth*(a1+ a2+ a3) + 9*seventh*nine_forty_cents*Area - Eu*invArea*(C12 + C13 + C23)/360.;
		invKBB = 1.0/KBB;

		// Vector Fh
		Fh[0] = twelfth*Area*(2*fe1 + fe2 + fe3);
		Fh[1] = twelfth*Area*(fe1 + 2*fe2 + fe3);
		Fh[2] = twelfth*Area*(fe1 + fe2 + 2*fe3);

		// Vector FB
		FB = three_twentieth*Area*(fe1 + fe2 + fe3);

		// Calculation of Ae and vector R
		for (i=0; i<3; i++){
			me[i][0] = Mhh_1_12  - invKBB*KhB[i]*MBh_3_20;
			me[i][1] = Mhh_1_12  - invKBB*KhB[i]*MBh_3_20;
			me[i][2] = Mhh_1_12  - invKBB*KhB[i]*MBh_3_20;
			me[i][i] += Mhh_1_12;	
			R[lm[E][i]] += Fh[i] - invKBB*KhB[i]*FB;
			for (j=0; j<3; j++){
				ke[i][j] = Khh[i][j] - invKBB*KhB[i]*KBh[j];
				Ae[i][j] = me[i][j] + alpha*dt*ke[i][j];
				R[lm[E][i]] -= me[i][j]*due[j] + ke[i][j]*ue[j];  		
			}
		}	
		R[NEQ] = 0;	
		
		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, E, Ae);
	

	}

	free(U);
	free(dU);

	return 0;
}

