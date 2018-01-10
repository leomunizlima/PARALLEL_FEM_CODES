#include "TranspEquation.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_M_K_R_SUPG(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int i, j, E, J1, J2, J3;
	int nel, NEQ;
	double kx, ky, Be_x, Be_y, gamma, TwoA, invArea, Area, norm, Be2_x, Be2_y, tau, h, h_shock, Eu, alpha, dt;
	double third=1.0/3.0, sixth = 1.0/6.0, twelfth= 1.0/12.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Beta[2], Ae[3][3], ke[3][3], me[3][3];
	double ue[3], ueb, due[3], dueb, fe1, fe2, fe3, feb, Rpg11, Rpg21, Rpg31, Cpg12, Cpg13, Cpg23 , Csc12, Csc13, Csc23;
	double D12, D13, D23, C11, C12, C13, R11, R12, XB, YB, *U, *dU;
	double *R = FemStructs->F;
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

	U = (double*) mycalloc("U of 'Build_M_F_DD_Transiente'", Parameters->nnodes, sizeof(double));
	dU = (double*) mycalloc("dU of 'Build_M_F_DD_Transiente'", Parameters->nnodes, sizeof(double));
	eval_U_dU(Parameters,FemStructs,FemFunctions,U,dU);

//===================================================================================================================================
/*	FILE *Out;
	char Filename[200];
	int a = (int) 100*Parameters->CurrentTime;
	sprintf(Filename,"RANK_%d_t%d.dat",Parameters->RANK,a);
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

		//Calculation of u and du in element e
		ue[0] = U[J1];
		ue[1] = U[J2];
		ue[2] = U[J3];
		due[0] = dU[J1];
		due[1] = dU[J2];
		due[2] = dU[J3];
		
		ueb = third * (ue[0] + ue[1] + ue[2]);
	   	dueb = third * (due[0] + due[1] + due[2]);

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

	   	
		/***************** local length element **********************/
		h_shock = FemFunctions->h_shock(Be_x, Be_y, ue[0], ue[1], ue[2], y23, y31, y12, x32, x13, x21, Area);
		/**************************************************************/

		/**************************** Shock capture parameter calculations***************************************/
		Eu = FemFunctions->ShockCapture(kx, ky, Be_x, Be_y, gamma, ue[0], ue[1], ue[2], ueb, dueb, feb, y23, y31,y12, x32,x13,x21,invArea, h_shock); 
		/********************************************************************************************************/

		// Shock capture correction matrix
		Csc12 = Eu*0.25*invArea*(y23*y31 + x32*x13);
	   	Csc13 = Eu*0.25*invArea*(y23*y12 + x32*x21);
	   	Csc23 = Eu*0.25*invArea*(y31*y12 + x13*x21);

	   	// Matrix Ae
	   	me[0][0] = sixth*Area  + tau*C11;
	   	me[0][1] = twelfth*Area + tau*C11;
	   	me[0][2] = twelfth*Area + tau*C11;
	    	me[1][0] = twelfth*Area + tau*C12;
	   	me[1][1] = sixth*Area  + tau*C12;
	    	me[1][2] = twelfth*Area + tau*C12;
	    	me[2][0] = twelfth*Area + tau*C13;
	    	me[2][1] = twelfth*Area + tau*C13;
	    	me[2][2] = sixth*Area  + tau*C13;
	
		ke[0][0] = -(D12 + D13) + C11 - (Cpg12 + Cpg13) + R11 + Rpg11 - (Csc12 + Csc13);
		ke[0][1] = 	    D12 + C12 + 	Cpg12 + R12 + Rpg11 + Csc12;
		ke[0][2] = 		D13 + C13 + 	Cpg13 + R12 + Rpg11 + Csc13;
		ke[1][0] = 		D12 + C11 + 	Cpg12 + R12 + Rpg21 + Csc12;
		ke[1][1] = -(D12 + D23) + C12 - (Cpg12 + Cpg23) + R11 + Rpg21 - (Csc12 + Csc23);
		ke[1][2] = 		D23 + C13 + Cpg23 + R12 + Rpg21 + Csc23;
		ke[2][0] = 		D13 + C11 + Cpg13 + R12 + Rpg31 + Csc13;
		ke[2][1] = 		D23 + C12 + Cpg23 + R12 + Rpg31 + Csc23;
		ke[2][2] = -(D13 + D23) + C13 - (Cpg13 + Cpg23) + R11 + Rpg31 - (Csc13 + Csc23);
		
		for (i=0; i<3; i++)
			for (j=0; j<3; j++)
				Ae[i][j] = me[i][j] + alpha*dt*ke[i][j];
		
		//Calculation of vector R 	
		R[lm[E][0]] += Area*twelfth*(2*fe1 +  fe2  +  fe3 ) + sixth*tau*(fe1+fe2+fe3)*(Be_x*y23 + Be_y*x32);
		R[lm[E][1]] += Area*twelfth*(  fe1 + 2*fe2 +  fe3 ) + sixth*tau*(fe1+fe2+fe3)*(Be_x*y31 + Be_y*x13);
		R[lm[E][2]] += Area*twelfth*(  fe1 +  fe2  + 2*fe3) + sixth*tau*(fe1+fe2+fe3)*(Be_x*y12 + Be_y*x21);

		int I;
		for (i=0; i<3; i++){
			I = lm[E][i];
			for (j=0; j<3; j++)
				R[I] -= me[i][j]*due[j] + ke[i][j]*ue[j];
			
		}
		R[NEQ] = 0;

/*		if (Parameters->RANK==0){
			if (E==1 || E==5 || E==6 || E==11){
				fprintf(Out,"E=%d x1:%lf x2:%lf x3:%lf y1:%lf y2:%lf y3:%lf\n",E, Node[J1].x,Node[J2].x,Node[J3].x,Node[J1].y,Node[J2].y,Node[J3].y);
				for (i=0;i<3;i++)
					for (j=0;j<3;j++)
						fprintf(Out,"me[%d][%d] = %.15lf ke[%d][%d]=%.15lf\n", i,j,me[i][j],i,j,ke[i][j]);

			}
		}

		if (Parameters->RANK==2){
			if (E==4 || E==5 || E==16 || E==17 || E==28 || E==29){*/
		/*		fprintf(Out,"E=%d x1:%lf x2:%lf x3:%lf y1:%lf y2:%lf y3:%lf\n",E, Node[J1].x,Node[J2].x,Node[J3].x,Node[J1].y,Node[J2].y,Node[J3].y);
				for (i=0;i<3;i++)
					for (j=0;j<3;j++)
						fprintf(Out,"me[%d][%d] = %.15lf ke[%d][%d]=%.15lf\n", i,j,me[i][j],i,j,ke[i][j]);
				for (i=0;i<3; i++)
					fprintf(Out,"R[%d]=%lf\n",lm[E][i],R[lm[E][i]]);

				fprintf(Out,"lm[%d][0]=%d lm[%d][1]=%d lm[%d][2]=%d \n",E,lm[E][0],E,lm[E][1],E,lm[E][2]);
	*/

//			}

//		}		
		// F contribution
	/*	double Fe1, Fe2, Fe3;
		Fe1 = Area*twelfth*(2*fe1 +  fe2  +  fe3 ) + sixth*tau*(fe1+fe2+fe3)*(Be_x*y23 + Be_y*x32);
		Fe2 = Area*twelfth*(  fe1 + 2*fe2 +  fe3 ) + sixth*tau*(fe1+fe2+fe3)*(Be_x*y31 + Be_y*x13);
		Fe3 = Area*twelfth*(  fe1 +  fe2  + 2*fe3) + sixth*tau*(fe1+fe2+fe3)*(Be_x*y12 + Be_y*x21);

		Fe1 -= me[0][0]*due[0] + ke[0][0]*ue[0] + me[0][1]*due[1] + ke[0][1]*ue[1] +me[0][2]*due[2] + ke[0][2]*ue[2];
		Fe2 -= me[1][0]*due[0] + ke[1][0]*ue[0] + me[1][1]*due[1] + ke[1][1]*ue[1] +me[1][2]*due[2] + ke[1][2]*ue[2];
		Fe3 -= me[2][0]*due[0] + ke[2][0]*ue[0] + me[2][1]*due[1] + ke[2][1]*ue[1] +me[2][2]*due[2] + ke[2][2]*ue[2];

		// Boundary conditions on vector F and Fill Matrix Diagonal 
		F_assembly(Node, J1, J2, J3, Fe1, Fe2, Fe3, R, ke, X, Y, FemFunctions);
*/
				
		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, E, Ae);

	}

//////////////////////////////////////////////////////////////
	
//	MPI_VectorUpdate(Parameters->RANK, Parameters->NPROC, Parameters, FemStructs, &R[Parameters->neq_bef]);
	

/*	int I;
	for (I=0; I<Parameters->nnodes; I++)
		fprintf(Out,"u[%d]=%lf du[%d]=%lf\n", I, U[I], I, dU[I]); 
	for (I=0; I<Parameters->nnodes; I++)
		if (Node[I].Type==1)
			fprintf(Out,"R[%d]=%lf\n", Node[I].invP_id, R[Node[I].invP_id+Parameters->neq_bef]); 

	fclose(Out);*/
	/////////////////////////////////////////////////////////////////////	

	
	free(U);
	free(dU);

	return 0;
}


