#include "TranspEquation.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_M_F_DD(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int i, j, E, J1, J2, J3;
	int nel, NEQ;
	double kx, ky, Be_x, Be_y, gamma, TwoA, invArea, Area, h_shock, Eu, alpha, dt;
	double w=0.5;
	double third=1.0/3.0, sixth = 1.0/6.0, seventh=1.0/7.0, twelfth= 1.0/12.0, three_twentieth=3./20.0, nine_twentieth=9./20.0, nine_forty_cents=9./40.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Beta[2], Me[3][3], XB, YB;
	double ue[3], ueb, due[3], dueb, fe1, fe2, fe3, feb;
	double **M2_out, **R2_out, *invN2_out, *U, *dU, *uB, *duB;
	double KDD12, KDD13, KDD23, D12, D13, D23, C12, C13, C23, R11, R12, a1, a2, a3, M1[3][3], N1[3], M2[3], N2, invN2, R1[3], R2, Mhh_1_12, MhB_3_20, MBh_3_20, MBB,
	Khh[3][3], KhB[3], KBh[3], KBB, Fh[3], FB, cbtil, cb;
	double *F = FemStructs->F;
	double *cbtilo = FemStructs->delta_old;
	int **lm = FemStructs->lm;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	nel = Parameters->nel;
	NEQ = Parameters->NEQ;
	kx = FemFunctions->Condutivity();
	ky = FemFunctions->Condutivity();
	gamma = FemFunctions->Reaction();
	M2_out = FemStructs->AuxBuild->M2;
	R2_out = FemStructs->AuxBuild->R2;
	invN2_out = FemStructs->AuxBuild->invN2;
	uB = FemStructs->uB;
	duB = FemStructs->duB;

	alpha = Parameters->Alpha_Build;
	dt = Parameters->DeltaT_Build;

	dzero(NEQ+1, F);
	setzeros(Parameters,MatrixData);

	U = (double*) mycalloc("U of 'Build_M_F_DD'", Parameters->nnodes, sizeof(double));
	dU = (double*) mycalloc("dU of 'Build_M_F_DD'", Parameters->nnodes, sizeof(double));
	eval_U_dU(Parameters,FemStructs,FemFunctions,U,dU);
//===================================================================================================================================
/*	FILE *Out;
	char Filename[200];
	int a = (int) 100*Parameters->CurrentTime;
	sprintf(Filename,"Data_%d_t%d.dat",Parameters->RANK,a);
	Out=fopen(Filename, "w");*/
//===================================================================================================================================

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

		//Matrix MBB
		MBB = (81.0/280.)*Area;
		
		//Matrix Khh
		//Termos auxiliares
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
	//	fprintf(Out,"cbtilo=%.20lf ",cbtilo[E]);
		cbtil = FemFunctions->ShockCapture(kx, ky, Be_x, Be_y, gamma, ue[0], ue[1], ue[2], ueb, dueb, feb, y23, y31,y12, x32,x13,x21,invArea, h_shock); 
	//	fprintf(Out,"cbtil=%.20lf\n",cbtil);
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
		Khh[1][2] = 		D23 + sixth*a3 + R12 + KDD23;
		Khh[2][0] = 		D13 + sixth*a1 + R12 + KDD13;
		Khh[2][1] = 		D23 + sixth*a2 + R12 + KDD23;
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

		// Vector Fh
		Fh[0] = twelfth*Area*(2*fe1 + fe2 + fe3);
		Fh[1] = twelfth*Area*(fe1 + 2*fe2 + fe3);
		Fh[2] = twelfth*Area*(fe1 + fe2 + 2*fe3);

		// Vector FB
		FB = three_twentieth*Area*(fe1 + fe2 + fe3);

		// Matrix N2 and invN2
		N2 = MBB + alpha*dt*KBB;
		invN2 = 1.0/N2;

		// Matrices M1, M2, N1, Me, and vectors R1, R2 calculations
		R2 = FB - MBB*duB[E] - KBB*uB[E];
		for (i=0; i<3; i++){
			M1[i][0] =  Mhh_1_12;
			M1[i][1] =  Mhh_1_12;
			M1[i][2] =  Mhh_1_12;
			M1[i][i] +=  Mhh_1_12;
			M2[i] = MBh_3_20 + alpha*dt*KBh[i];
			N1[i] = MhB_3_20 + alpha*dt*KhB[i];
			R1[i] = Fh[i] - MhB_3_20*duB[E] - KhB[i]*uB[E];
			R1[i] -= Mhh_1_12*due[i];
			for (j=0; j<3; j++){
				M1[i][j] +=  alpha*dt*Khh[i][j];
				R1[i] -= Mhh_1_12*due[j] + Khh[i][j]*ue[j];
			}
			R2 -= MBh_3_20*due[i]+ KBh[i]*ue[i];
		}

		for (i=0; i<3; i++){
			for (j=0; j<3; j++)
				Me[i][j] = M1[i][j] - invN2*N1[i]*M2[j];
			F[lm[E][i]] += R1[i] - invN2*N1[i]*R2;
		}
		F[NEQ] = 0;
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
/*		fprintf(Out,"E=%d x1:%lf x2:%lf x3:%lf y1:%lf y2:%lf y3:%lf\n",E, Node[J1].x,Node[J2].x,Node[J3].x,Node[J1].y,Node[J2].y,Node[J3].y);
		fprintf(Out,"invN2=%.20lf\n",invN2);
		fprintf(Out,"R2=%.20lf\n",R2);
		fprintf(Out,"ue[0]=%.lf ue[1]=%.20lf ue[2]=%.20lf\n",ue[0],ue[1],ue[2]);
		fprintf(Out,"due[0]=%.lf due[1]=%.20lf due[2]=%.20lf\n",due[0],due[1],due[2]);
		fprintf(Out,"uB=%.20lf duB=%.20lf\n",uB[E], duB[E]);
		for (i=0;i<3;i++)
			fprintf(Out,"R1[%d]=%.20lf\n",i,R1[i]);
		for (i=0;i<3;i++)
			fprintf(Out,"N1[%d]=%.20lf\n",i,N1[i]);
		for (i=0;i<3;i++)
			for (j=0;j<3;j++)
				fprintf(Out,"Me[%d][%d] = %.20lf\n", i,j,Me[i][j]);
		for (i=0;i<3; i++)
			fprintf(Out,"F[%d]=%.20lf\n",lm[E][i],F[lm[E][i]]);
		
*/	
//-------------------------------------------------------------------------------------------------------------------------------------------------------

		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, E, Me);
	

		// Structs for DduB calculation
		invN2_out[E] = invN2; 
		R2_out[E][0] = R2;
		M2_out[E][0] = M2[0];
		M2_out[E][1] = M2[1];
		M2_out[E][2] = M2[2];
				
	}
	
//-------------------------------------------------------------------------------------------------------------------------------------------------------
/*	int I;
	for (I=0; I<Parameters->nnodes; I++)
		if (Node[I].Type==1)
			fprintf(Out,"u[%d]=%.20lf du[%d]=%.20lf\n", Node[I].invP_id, FemStructs->u[Node[I].invP_id], Node[I].invP_id, FemStructs->du[Node[I].invP_id]);
	for (I=0; I<Parameters->nnodes; I++)
		if (Node[I].Type==1)
			fprintf(Out,"F[%d]=%.20lf\n", Node[I].invP_id, F[Node[I].invP_id+Parameters->neq_bef]); 

	fclose(Out);*/
//-------------------------------------------------------------------------------------------------------------------------------------------------------

	free(U);
	free(dU);

	return 0;
}


