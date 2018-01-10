#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

static int varglobal = 0;

int Build_K_F_VMS_DCDD(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int i, NEQ, nel;//, op;
	int E, J1, J2, J3;
	//double Uref, Lref;
	double Re, x, y;
	double TwoA, invArea, Area, tau, h, visc, rho, unorm, tau_l, normgradU; //, invrho
	double third=1.0/3.0, sixth = 1.0/6.0,  ninefortieth = 9.0 / 40.0; //, twelfth = 1.0/12.0, ninetwentieth = 9.0/20.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Ke[9][9], Fe[9], ux1, ux2, ux3, uy1, uy2, uy3, ux, uy, duxdx, duxdy, duydx, duydy;
	double fx1, fx2, fx3, fy1, fy2, fy3, fxB, fyB, f1, f2, fdelta1, fdelta2, fdelta3, fdelta4, fdelta5, fdelta6, fphi1, fphi2, fphi3;
	double C1, C2, C3;
	double Nv[6], Ndeltav[6], Kv[6], Ksv[6], Gv[6], Gdeltav[6], GTv[3], Nphiv[3], Gphiv[3], Res[9];
	//double M11, M13, Mdelta11, Mdelta31, Mdelta51;
	double K11, K12, K13, K14, K15, K16, K22, K23, K24, K25, K26, K33, K34, K35, K36, K44, K45, K46, K55, K56, K66; 
	double Ks11, Ks12, Ks13, Ks14, Ks15, Ks16, Ks21, Ks22, Ks23, Ks24, Ks25, Ks26, Ks31, Ks32, Ks33, Ks34, Ks35, Ks36;
	double Ks41, Ks42, Ks43, Ks44, Ks45, Ks46, Ks51, Ks52, Ks53, Ks54, Ks55, Ks56, Ks61, Ks62, Ks63, Ks64, Ks65, Ks66; 
	double GT11, GT12, GT13, GT14, GT15, GT16, N11, N13, N15;
	double Ndelta11, Ndelta12, Ndelta13, Ndelta14, Ndelta15, Ndelta16, Ndelta21, Ndelta22, Ndelta23, Ndelta24, Ndelta25, Ndelta26;
	double Ndelta31, Ndelta32, Ndelta33, Ndelta34, Ndelta35, Ndelta36, Ndelta41, Ndelta42, Ndelta43, Ndelta44, Ndelta45, Ndelta46;
	double Ndelta51, Ndelta52, Ndelta53, Ndelta54, Ndelta55, Ndelta56, Ndelta61, Ndelta62, Ndelta63, Ndelta64, Ndelta65, Ndelta66;
	double Gdelta11, Gdelta12, Gdelta13, Gdelta21, Gdelta22, Gdelta23, Gdelta31, Gdelta32, Gdelta33, Gdelta41, Gdelta42, Gdelta43;
	double Gdelta51, Gdelta52, Gdelta53, Gdelta61, Gdelta62, Gdelta63;
	//double Mphi11, Mphi12, Mphi21, Mphi22, Mphi31, Mphi32;
	double Nphi11, Nphi12, Nphi13, Nphi14, Nphi15, Nphi16, Nphi21, Nphi22, Nphi23, Nphi24, Nphi25, Nphi26, Nphi31, Nphi32, Nphi33;
	double Nphi34, Nphi35, Nphi36;	
	double Gphi11, Gphi12, Gphi13, Gphi21, Gphi22, Gphi23, Gphi31, Gphi32, Gphi33;
	//double Np11, Np12, Np21, Np22, Np31, Np32, Np41, Np42;
	//double Npdelta11, Npdelta12, Npdelta21, Npdelta22, Npdelta31, Npdelta32, Npdelta41, Npdelta42, Npdelta51, Npdelta52, Npdelta61, Npdelta62; 
	//double Nppdelta11, Nppdelta12, Nppdelta21, Nppdelta22, Nppdelta31, Nppdelta32, Nppdelta41, Nppdelta42, Nppdelta51, Nppdelta52, Nppdelta61, Nppdelta62;
	//double Npphi11, Npphi12, Npphi21, Npphi22, Npphi31, Npphi32;
	double p1, p2, p3, p, duxduy;
	
	double *F = FemStructs->F;
	int **lm = FemStructs->lm;
	
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	nel = Parameters->nel;
	NEQ = Parameters->NEQ;
	
	dzero(NEQ+1, F);
	setzeros(Parameters,MatrixData);

	double *U = (double*) mycalloc("U of 'Build_K_F_VMS_DCDD'", 3*Parameters->nnodes, sizeof(double));
	eval_U(Parameters,FemStructs,FemFunctions,U);

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

		TwoA =  fabs(x21*y31 - x13*y12);
		Area = 0.5*TwoA;
		invArea = 1.0/Area;		
				
		//x Velocity  
		ux1 = U[3*J1];
		ux2 = U[3*J2];
		ux3 = U[3*J3];

		//y Velocity
		uy1 = U[3*J1+1];
		uy2 = U[3*J2+1];
		uy3 = U[3*J3+1];
		
		//pression 
		p1 = U[3*J1+2];
		p2 = U[3*J2+2];
		p3 = U[3*J3+2];
		
		
		//*****************u at barycentre of the element***************
		ux = third*(ux1 + ux2 + ux3);
		uy = third*(uy1 + uy2 + uy3);
		
		//*****************************Gradient of u********************
		duxdx = 0.5*invArea*( ux1*y23 + ux2*y31 + ux3*y12 ); 
		duxdy = 0.5*invArea*( ux1*x32 + ux2*x13 + ux3*x21 );
		duydx = 0.5*invArea*( uy1*y23 + uy2*y31 + uy3*y12 );
		duydy = 0.5*invArea*( uy1*x32 + uy2*x13 + uy3*x21 );
			
		duxduy = ux1*y23 + uy1*x32 + ux2*y31 + uy2*x13 + ux3*y12 + uy3*x21;


		//*****************p at barycentre of the element***************
		p = third*( p1 + p2 + p3 );

		//**************************Gradient of p**********************
		//dpdx = 0.5*invArea*( p1*y23 + p2*y31 + p3*y12 );
		//dpdy = 0.5*invArea*( p1*x32 + p2*x13 + p3*x21 );

		//*****************Coefficients C_i******************************
		C1 = 0.5*invArea*( ux*y23 + uy*x32 );
		C2 = 0.5*invArea*( ux*y31 + uy*x13 );
		C3 = 0.5*invArea*( ux*y12 + uy*x21 );

		//************* Dados de entrada *********
		Re = Parameters->ReynoldsNumber;		
		h = sqrt( 4*Area/PI );		
		visc =1./Re;
		rho = 1.;		
		//invrho = 1.0/rho;
		//nu = visc/rho;
		//h = sqrt( 4*Area/PI );		
		unorm = sqrt( ux*ux + uy*uy );
		
		//auxtau = ( 2*unorm/h )*( 2*unorm/h ) + 9*( 4*nu/(h*h) )*( 4*nu/(h*h) );
		//auxtau = sqrt( auxtau);
		//tau = 1.0/auxtau;
		//tau=tau_SUPG=tau_PSPG
		
		//*************Calculation of tau_LSIC stabilization*********
		tau_l = 0.5*h*unorm ;
		
	
		//************Matrices of the Galerkin formulation*************

		//*************************Matrix M****************************
		//M11 = rho*2*Area*twelfth;
		//M13 = rho*Area*twelfth; 
		
		//*************************Matrix N****************************
		N11 = rho*third*Area*C1;		
		N13 = rho*third*Area*C2;		
		N15 = rho*third*Area*C3;		

		//************************Matrix K******************************
		K11 = visc*0.25*invArea*( 2*y23*y23 + x32*x32 );
		K12 = visc*0.25*invArea*( x32*y23 );
		K13 = visc*0.25*invArea*( 2*y23*y31 + x32*x13 );
		K14 = visc*0.25*invArea*( x32*y31 );
		K15 = visc*0.25*invArea*( 2*y23*y12 + x32*x21 );
		K16 = visc*0.25*invArea*( x32*y12 );

		K22 = visc*0.25*invArea*( y23*y23 + 2*x32*x32 );
		K23 = visc*0.25*invArea*( y23*x13 );
		K24 = visc*0.25*invArea*( y23*y31 + 2*x32*x13 );
		K25 = visc*0.25*invArea*( y23*x21 );
		K26 = visc*0.25*invArea*( y23*y12 + 2*x32*x21 );

		K33 = visc*0.25*invArea*( 2*y31*y31 + x13*x13 );
		K34 = visc*0.25*invArea*( x13*y31 );
		K35 = visc*0.25*invArea*( 2*y31*y12 + x13*x21 );
		K36 = visc*0.25*invArea*( x13*y12 );

		K44 = visc*0.25*invArea*( y31*y31 + 2*x13*x13 );
		K45 = visc*0.25*invArea*( y31*x21 );
		K46 = visc*0.25*invArea*( y31*y12 + 2*x13*x21 );

		K55 = visc*0.25*invArea*( 2*y12*y12 + x21*x21 );		
		K56 = visc*0.25*invArea*( x21*y12 );

		K66 = visc*0.25*invArea*( y12*y12 + 2*x21*x21);

		//************************Matrix G^T****************************
		GT11 = sixth*y23;
		GT12 = sixth*x32;
		GT13 = sixth*y31;
		GT14 = sixth*x13;
		GT15 = sixth*y12;
		GT16 = sixth*x21;
		
		//***************Matrices of the LSIC formulation****************
		
		//***************Matrices K_s ***********************************
		
		Ks11 = tau_l*rho*0.25*invArea*y23*y23;
		Ks12 = tau_l*rho*0.25*invArea*x32*y23;
		Ks13 = tau_l*rho*0.25*invArea*y31*y23;
		Ks14 = tau_l*rho*0.25*invArea*x13*y23;
		Ks15 = tau_l*rho*0.25*invArea*y12*y23;
		Ks16 = tau_l*rho*0.25*invArea*x21*y23;

		Ks21 = tau_l*rho*0.25*invArea*y23*x32;
		Ks22 = tau_l*rho*0.25*invArea*x32*x32;
		Ks23 = tau_l*rho*0.25*invArea*y31*x32;
		Ks24 = tau_l*rho*0.25*invArea*x13*x32;
		Ks25 = tau_l*rho*0.25*invArea*y12*x32;
		Ks26 = tau_l*rho*0.25*invArea*x21*x32;

		Ks31 = tau_l*rho*0.25*invArea*y23*y31;
		Ks32 = tau_l*rho*0.25*invArea*x32*y31;
		Ks33 = tau_l*rho*0.25*invArea*y31*y31;
		Ks34 = tau_l*rho*0.25*invArea*x13*y31;
		Ks35 = tau_l*rho*0.25*invArea*y12*y31;
		Ks36 = tau_l*rho*0.25*invArea*x21*y31;

		Ks41 = tau_l*rho*0.25*invArea*y23*x13;
		Ks42 = tau_l*rho*0.25*invArea*x32*x13;
		Ks43 = tau_l*rho*0.25*invArea*y31*x13;
		Ks44 = tau_l*rho*0.25*invArea*x13*x13;
		Ks45 = tau_l*rho*0.25*invArea*y12*x13;
		Ks46 = tau_l*rho*0.25*invArea*x21*x13;
		
		Ks51 = tau_l*rho*0.25*invArea*y23*y12;
		Ks52 = tau_l*rho*0.25*invArea*x32*y12;
		Ks53 = tau_l*rho*0.25*invArea*y31*y12;
		Ks54 = tau_l*rho*0.25*invArea*x13*y12;
		Ks55 = tau_l*rho*0.25*invArea*y12*y12;
		Ks56 = tau_l*rho*0.25*invArea*x21*y12;

		Ks61 = tau_l*rho*0.25*invArea*y23*x21;
		Ks62 = tau_l*rho*0.25*invArea*x32*x21;
		Ks63 = tau_l*rho*0.25*invArea*y31*x21;
		Ks64 = tau_l*rho*0.25*invArea*x13*x21;
		Ks65 = tau_l*rho*0.25*invArea*y12*x21;
		Ks66 = tau_l*rho*0.25*invArea*x21*x21;
		

		
		//*********************Incremental matrices***********************
		
		
		//****************************Font***********************************
/*		fx1 = (*Font)(X[0],Y[0],1);*/
/*		fx2 = (*Font)(X[1],Y[1],1);*/
/*		fx3 = (*Font)(X[2],Y[2],1);*/
/*		*/
/*		fy1 = (*Font)(X[0],Y[0],2);*/
/*		fy2 = (*Font)(X[1],Y[1],2);*/
/*		fy3 = (*Font)(X[2],Y[2],2);*/

		fx1 = 0.0;
		fx2 = 0.0;
		fx3 = 0.0;
		
		fy1 = 0.0;
		fy2 = 0.0;
		fy3 = 0.0;


		//*********************Fonte Sol Exata Conhecida********************
		x = X[0];
		y = Y[0];
		
		fx1 = 2*x + pow(x,4)*(1.2 - 2.4*y) + pow(x,3)*(-2.4 + 4.8*y) + 
 y*(-0.4 + 1.2*y - 0.8*pow(y,2)) + 8*pow((-1 + x),3)*pow(x,3)*(-1 + 2*x)*pow(y,2)*pow((1 - 3*y + 2*pow(y,2)),2) + x*y*(2.4 - 7.2*y + 4.8*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,3)*(1 - 3*x + 2*pow(x,2))*pow((-1 + y),2)*pow(y,2)*(1 - 6*y + 6*pow(y,2)) + pow(x,2)*(1.2 - 4.8*y + 7.2*pow(y,2) - 4.8*pow(y,3));
		
		fy1 =	-2*y - 1.2*pow((1 - y),2)*pow(y,2) + 
 8*pow(x,2)*pow((1 - 3*x + 2*pow(x,2)),2)*pow((-1 + y),3)*pow(y,3)*(-1 + 2*y) + 
 pow(x,2)*(-1.2 + 7.2*y - 7.2*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,2)*(1 - 6*x + 6*pow(x,2))*pow((-1 + y),2)*pow(y,3)*(1 - 3*y + 2*pow(y,2)) + pow(x,3)*(0.8 - 4.8*y + 4.8*pow(y,2)) + x*(0.4 - 2.4*y + 4.8*pow(y,2) - 4.8*pow(y,3) + 2.4*pow(y,4));	

		x = X[1];
		y = Y[1];
		
		fx2 = 2*x + pow(x,4)*(1.2 - 2.4*y) + pow(x,3)*(-2.4 + 4.8*y) + 
 y*(-0.4 + 1.2*y - 0.8*pow(y,2)) + 8*pow((-1 + x),3)*pow(x,3)*(-1 + 2*x)*pow(y,2)*pow((1 - 3*y + 2*pow(y,2)),2) + x*y*(2.4 - 7.2*y + 4.8*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,3)*(1 - 3*x + 2*pow(x,2))*pow((-1 + y),2)*pow(y,2)*(1 - 6*y + 6*pow(y,2)) + pow(x,2)*(1.2 - 4.8*y + 7.2*pow(y,2) - 4.8*pow(y,3));
		
		fy2 =	-2*y - 1.2*pow((1 - y),2)*pow(y,2) + 
 8*pow(x,2)*pow((1 - 3*x + 2*pow(x,2)),2)*pow((-1 + y),3)*pow(y,3)*(-1 + 2*y) + 
 pow(x,2)*(-1.2 + 7.2*y - 7.2*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,2)*(1 - 6*x + 6*pow(x,2))*pow((-1 + y),2)*pow(y,3)*(1 - 3*y + 2*pow(y,2)) + pow(x,3)*(0.8 - 4.8*y + 4.8*pow(y,2)) + x*(0.4 - 2.4*y + 4.8*pow(y,2) - 4.8*pow(y,3) + 2.4*pow(y,4));
		
		x = X[2];
		y = Y[2];
		
		fx3 = 2*x + pow(x,4)*(1.2 - 2.4*y) + pow(x,3)*(-2.4 + 4.8*y) + 
 y*(-0.4 + 1.2*y - 0.8*pow(y,2)) + 8*pow((-1 + x),3)*pow(x,3)*(-1 + 2*x)*pow(y,2)*pow((1 - 3*y + 2*pow(y,2)),2) + x*y*(2.4 - 7.2*y + 4.8*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,3)*(1 - 3*x + 2*pow(x,2))*pow((-1 + y),2)*pow(y,2)*(1 - 6*y + 6*pow(y,2)) + pow(x,2)*(1.2 - 4.8*y + 7.2*pow(y,2) - 4.8*pow(y,3));
		
		fy3 =	-2*y - 1.2*pow((1 - y),2)*pow(y,2) + 
 8*pow(x,2)*pow((1 - 3*x + 2*pow(x,2)),2)*pow((-1 + y),3)*pow(y,3)*(-1 + 2*y) + 
 pow(x,2)*(-1.2 + 7.2*y - 7.2*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,2)*(1 - 6*x + 6*pow(x,2))*pow((-1 + y),2)*pow(y,3)*(1 - 3*y + 2*pow(y,2)) + pow(x,3)*(0.8 - 4.8*y + 4.8*pow(y,2)) + x*(0.4 - 2.4*y + 4.8*pow(y,2) - 4.8*pow(y,3) + 2.4*pow(y,4));

		fxB = third*( fx1 + fx2 + fx3 );
		fyB = third*( fy1 + fy2 + fy3 );
		
		//***************Font of the Galerkin formulation********************
		f1 = rho*third*Area*fxB;
		f2 = rho*third*Area*fyB;

		f1 = f2 = fxB = fyB = 0.0;
		
		//*****************Font of the SUPG formulation**********************
		fdelta1 = rho*tau*Area*C1*fxB;
		fdelta2 = rho*tau*Area*C1*fyB;
		fdelta3 = rho*tau*Area*C2*fxB;
		fdelta4 = rho*tau*Area*C2*fyB;
		fdelta5 = rho*tau*Area*C3*fxB;
		fdelta6 = rho*tau*Area*C3*fyB;

		//*****************Font of the PSPG formulation**********************
		fphi1 = tau*Area*( y23*fxB + x32*fyB );
		fphi2 = tau*Area*( y31*fxB + x13*fyB );
		fphi3 = tau*Area*( y12*fxB + x21*fyB );

		//*******************************vector Fe***************************
		Fe[0] = f1 + fdelta1;
		Fe[1] = f2 + fdelta2;
		Fe[2] = fphi1;
		Fe[3] = f1 + fdelta3;
		Fe[4] = f2 + fdelta4;
		Fe[5] = fphi2;
		Fe[6] = f1 + fdelta5;
		Fe[7] = f2 + fdelta6;
		Fe[8] = fphi3;

		//**********************Calculation of the residue*******************		
		Nv[0] = rho*third*Area*( duxdx*ux + duxdy*uy );
		Nv[1] = rho*third*Area*( duydx*ux + duydy*uy );
		Nv[2] = rho*third*Area*( duxdx*ux + duxdy*uy );
		Nv[3] = rho*third*Area*( duydx*ux + duydy*uy );
		Nv[4] = rho*third*Area*( duxdx*ux + duxdy*uy );
		Nv[5] = rho*third*Area*( duydx*ux + duydy*uy );

		Kv[0] = 0.5*visc*( 2*y23*duxdx + x32*( duxdy + duydx ) );
		Kv[1] = 0.5*visc*( 2*x32*duydy + y23*( duxdy + duydx ) );
		Kv[2] = 0.5*visc*( 2*y31*duxdx + x13*( duxdy + duydx ) );
		Kv[3] = 0.5*visc*( 2*x13*duydy + y31*( duxdy + duydx ) );
		Kv[4] = 0.5*visc*( 2*y12*duxdx + x21*( duxdy + duydx ) );
		Kv[5] = 0.5*visc*( 2*x21*duydy + y12*( duxdy + duydx ) );
		
		Ksv[0] = tau_l*rho*0.25*invArea*y23*duxduy;
		Ksv[1] = tau_l*rho*0.25*invArea*x32*duxduy;
		Ksv[2] = tau_l*rho*0.25*invArea*y31*duxduy;
		Ksv[3] = tau_l*rho*0.25*invArea*x13*duxduy;
		Ksv[4] = tau_l*rho*0.25*invArea*y12*duxduy;
		Ksv[5] = tau_l*rho*0.25*invArea*x21*duxduy;

		Gv[0] = 0.5*p*y23; 
		Gv[1] = 0.5*p*x32;
		Gv[2] = 0.5*p*y31;
		Gv[3] = 0.5*p*x13;
		Gv[4] = 0.5*p*y12;
		Gv[5] = 0.5*p*x21;

		GTv[0] = third*Area*( duxdx + duydy);
		GTv[1] = third*Area*( duxdx + duydy);
		GTv[2] = third*Area*( duxdx + duydy);

		//*****************Residuo para o calculo do coeficiente da difusao artificial 
		
		Res[0] = Fe[0] - (Nv[0] + Kv[0] - Gv[0]);
		Res[1] = Fe[1] - (Nv[1] + Kv[1] - Gv[1]);
		Res[2] = Fe[2] - GTv[0] ;
		Res[3] = Fe[3] - (Nv[2] + Kv[2] - Gv[2] );
		Res[4] = Fe[4] - (Nv[3] + Kv[3] - Gv[3] );
		Res[5] = Fe[5] - GTv[1] ;
		Res[6] = Fe[6] - (Nv[4] + Kv[4] - Gv[4]);
		Res[7] = Fe[7] - (Nv[5] + Kv[5] - Gv[5]);
		Res[8] = Fe[8] - GTv[2];

		//*****************Norma euclidiana do residuo ||R(U,P)||		
	//	normres = sqrt(Res[0]*Res[0] + Res[1]*Res[1] + Res[2]*Res[2] + Res[3]*Res[3] + Res[4]*Res[4] + Res[5]*Res[5] + Res[6]*Res[6] + Res[7]*Res[7] + Res[8]*Res[8]);

		//*****************Norma do grad de U
		normgradU = sqrt(duxdx*duxdx + duxdy*duxdy + duydx*duydx + duydy*duydy);
		
		//******Analogous matrices of the SUPG formulation***************

		//*******************Matrix Mdelta=tau*N^T***********************
		//Mdelta11 = tau*N11;
		//Mdelta31 = tau*N13;
		//Mdelta51 = tau*N15;

		//*********Matrix Ndelta = NhB(KBB + KBBdelta)^{-1}NBh**********
	
			//*************** Matrix NBh ***************************
		double NBh11, NBh13, NBh15;
		NBh11 = rho*ninefortieth*(ux*y23 + uy*x32);
		NBh13 = rho*ninefortieth*(ux*y31 + uy*x13);
		NBh15 = rho*ninefortieth*(ux*y12 + uy*x21);
		

			//************** Matrix NhB = - NBh^T *******************		
		double NhB11, NhB31, NhB51;
		NhB11 = -NBh11;
		NhB31 = -NBh13;
		NhB51 = -NBh15;

			//****** Matrix InvKBBd = (KBB + KBBdelta)^{-1}**********
		double A, B, C, D, InvKBBd11, InvKBBd12, InvKBBd22, delta;
		A = y23*y23 + y31*y31 + y12*y12 + y23*y31 + y23*y12 + y31*y12;
		B = x32*x32 + x13*x13 + x21*x21 + x32*x13 + x32*x21 + x13*x21;
		C = y23*x32 + y31*x13 + y12*x21;

			//******* Parametro de difusao artificial := delta ******		
		if(normgradU>0.0){		
			delta = 0.0;//0.5*h*normres/(2.0*normgradU);
			//printf( "Delta NAO nulo = %15.14f \n", delta);
		}else{
			delta = 0.0;
			//printf( "Delta nulo\n");
		}
	
		//******** D eh o determinande da matriz (KBB + KBBdelta) 
		D = ((2*visc+delta)*(2*visc+delta) + (visc+delta)*(visc+delta))*A*B + (2*visc+delta)*(visc+delta)*(A*A + B*B) - visc*visc*C*C/4.0;
		
		InvKBBd11 = 40.0*Area/(81.0*D)*((2*visc+delta)*B + (visc+delta)*A);
		InvKBBd12 = - 40.0*Area/(81.0*D)*(visc*C/2.0);
		InvKBBd22 = 40.0*Area/(81.0*D)*((2*visc+delta)*A + (visc+delta)*B);

		Ndelta11 = NhB11*InvKBBd11*NBh11;
		Ndelta12 = NhB11*InvKBBd12*NBh11;
		Ndelta13 = NhB11*InvKBBd11*NBh13;
		Ndelta14 = NhB11*InvKBBd12*NBh13;
		Ndelta15 = NhB11*InvKBBd11*NBh15;
		Ndelta16 = NhB11*InvKBBd12*NBh15;
		
		Ndelta21 = Ndelta12;
		Ndelta22 = NhB11*InvKBBd22*NBh11;
		Ndelta23 = Ndelta14; 		
		Ndelta24 = NhB11*InvKBBd22*NBh13;
		Ndelta25 = Ndelta16; 		
		Ndelta26 = NhB11*InvKBBd22*NBh15;

		Ndelta31 = NhB31*InvKBBd11*NBh11;
		Ndelta32 = NhB31*InvKBBd12*NBh11;
		Ndelta33 = NhB31*InvKBBd11*NBh13;
		Ndelta34 = NhB31*InvKBBd12*NBh13;
		Ndelta35 = NhB31*InvKBBd11*NBh15;
		Ndelta36 = NhB31*InvKBBd12*NBh15;
		
		Ndelta41 = Ndelta32;
		Ndelta42 = NhB31*InvKBBd22*NBh11;
		Ndelta43 = Ndelta34;
		Ndelta44 = NhB31*InvKBBd22*NBh13;
		Ndelta45 = Ndelta36;		
		Ndelta46 = NhB31*InvKBBd22*NBh15;

		Ndelta51 = NhB51*InvKBBd11*NBh11;
		Ndelta52 = NhB51*InvKBBd12*NBh11;
		Ndelta53 = NhB51*InvKBBd11*NBh13;
		Ndelta54 = NhB51*InvKBBd12*NBh13;
		Ndelta55 = NhB51*InvKBBd11*NBh15;
		Ndelta56 = NhB51*InvKBBd12*NBh15;
	
		Ndelta61 = Ndelta52;
		Ndelta62 = NhB51*InvKBBd22*NBh11;
		Ndelta63 = Ndelta54;
		Ndelta64 = NhB51*InvKBBd22*NBh13;
		Ndelta65 = Ndelta56;
		Ndelta66 = NhB51*InvKBBd22*NBh15;


		//*********Matrix Gdelta = NhB(KBB + KBBdelta)^{-1}GBh**********
			//************** Matrix GBh ****************************
		double GBh11, GBh12, GBh13, GBh21, GBh22, GBh23;		
		GBh11 = - ninefortieth*y23;		
		GBh12 = - ninefortieth*y31;		
		GBh13 = - ninefortieth*y12;		
		GBh21 = - ninefortieth*x32;		
		GBh22 = - ninefortieth*x13;		
		GBh23 = - ninefortieth*x21;		
			//*************Matrix Gdelta****************************
		Gdelta11 = NhB11*InvKBBd11*GBh11 + NhB11*InvKBBd12*GBh21;
		Gdelta12 = NhB11*InvKBBd11*GBh12 + NhB11*InvKBBd12*GBh22;
		Gdelta13 = NhB11*InvKBBd11*GBh13 + NhB11*InvKBBd12*GBh23;
		
		Gdelta21 = NhB11*InvKBBd12*GBh11 + NhB11*InvKBBd22*GBh21;
		Gdelta22 = NhB11*InvKBBd12*GBh12 + NhB11*InvKBBd22*GBh22;
		Gdelta23 = NhB11*InvKBBd12*GBh13 + NhB11*InvKBBd22*GBh23;
			
		Gdelta31 = NhB31*InvKBBd11*GBh11 + NhB31*InvKBBd12*GBh21;
		Gdelta32 = NhB31*InvKBBd11*GBh12 + NhB31*InvKBBd12*GBh22;
		Gdelta33 = NhB31*InvKBBd11*GBh13 + NhB31*InvKBBd12*GBh23;
		
		Gdelta41 = NhB31*InvKBBd12*GBh11 + NhB31*InvKBBd22*GBh21;
		Gdelta42 = NhB31*InvKBBd12*GBh12 + NhB31*InvKBBd22*GBh22;
		Gdelta43 = NhB31*InvKBBd12*GBh13 + NhB31*InvKBBd22*GBh23;
		
		Gdelta51 = NhB51*InvKBBd11*GBh11 + NhB51*InvKBBd12*GBh21;
		Gdelta52 = NhB51*InvKBBd11*GBh12 + NhB51*InvKBBd12*GBh22;
		Gdelta53 = NhB51*InvKBBd11*GBh13 + NhB51*InvKBBd12*GBh23;
		
		Gdelta61 = NhB51*InvKBBd12*GBh11 + NhB51*InvKBBd22*GBh21;
		Gdelta62 = NhB51*InvKBBd12*GBh12 + NhB51*InvKBBd22*GBh22;
		Gdelta63 = NhB51*InvKBBd12*GBh13 + NhB51*InvKBBd22*GBh23;
		
		//*******Analogous matrices of the PSPG formulation**************

		//**********************Matrix Mphi******************************
		//Mphi11 = tau*sixth*y23;
		//Mphi12 = tau*sixth*x32;
	
		//Mphi21 = tau*sixth*y31;
		//Mphi22 = tau*sixth*x13;
	
		//Mphi31 = tau*sixth*y12;
		//Mphi32 = tau*sixth*x21;

		//** Matrix Nphi = GhB(KBB+KBBdelta)^{-1}NBh = - Gdelta^T*************
		Nphi11 = - Gdelta11;		
		Nphi12 = - Gdelta21;		
		Nphi13 = - Gdelta31;		
		Nphi14 = - Gdelta41;		
		Nphi15 = - Gdelta51;		
		Nphi16 = - Gdelta61;		
		
		Nphi21 = - Gdelta12;		
		Nphi22 = - Gdelta22;		
		Nphi23 = - Gdelta32;		
		Nphi24 = - Gdelta42;		
		Nphi25 = - Gdelta52;		
		Nphi26 = - Gdelta62;		
		
		Nphi31 = - Gdelta13;		
		Nphi32 = - Gdelta23;		
		Nphi33 = - Gdelta33;		
		Nphi34 = - Gdelta43;		
		Nphi35 = - Gdelta53;		
		Nphi36 = - Gdelta63;		
		
		//** Matrix Gphi = GhB(KBB+KBBdelta)^{-1}GBh, with GhB = GBh^T*************
 	
		Gphi11 = GBh11*(GBh11*InvKBBd11 + GBh21*InvKBBd12) + GBh21*(GBh11*InvKBBd12 + GBh21*InvKBBd22);
		Gphi12 = GBh12*(GBh11*InvKBBd11 + GBh21*InvKBBd12) + GBh22*(GBh11*InvKBBd12 + GBh21*InvKBBd22);
		Gphi13 = GBh13*(GBh11*InvKBBd11 + GBh21*InvKBBd12) + GBh23*(GBh11*InvKBBd12 + GBh21*InvKBBd22);
		
		Gphi21 = GBh11*(GBh12*InvKBBd11 + GBh22*InvKBBd12) + GBh21*(GBh12*InvKBBd12 + GBh22*InvKBBd22);
		Gphi22 = GBh12*(GBh12*InvKBBd11 + GBh22*InvKBBd12) + GBh22*(GBh12*InvKBBd12 + GBh22*InvKBBd22);
		Gphi23 = GBh13*(GBh12*InvKBBd11 + GBh22*InvKBBd12) + GBh23*(GBh12*InvKBBd12 + GBh22*InvKBBd22);
		
		Gphi31 = GBh11*(GBh13*InvKBBd11 + GBh23*InvKBBd12) + GBh21*(GBh13*InvKBBd12 + GBh23*InvKBBd22);
		Gphi32 = GBh12*(GBh13*InvKBBd11 + GBh23*InvKBBd12) + GBh22*(GBh13*InvKBBd12 + GBh23*InvKBBd22);
		Gphi33 = GBh13*(GBh13*InvKBBd11 + GBh23*InvKBBd12) + GBh23*(GBh13*InvKBBd12 + GBh23*InvKBBd22);
		
		//*************Calculation of the residue continue*******************

		Ndeltav[0] = Ndelta11*ux1 +  Ndelta12*uy1 +  Ndelta13*ux2 +  Ndelta14*uy2 +  Ndelta15*ux3 +  Ndelta16*uy3;
		Ndeltav[1] = Ndelta21*ux1 +  Ndelta22*uy1 +  Ndelta23*ux2 +  Ndelta24*uy2 +  Ndelta25*ux3 +  Ndelta26*uy3;
		Ndeltav[2] = Ndelta31*ux1 +  Ndelta32*uy1 +  Ndelta33*ux2 +  Ndelta34*uy2 +  Ndelta35*ux3 +  Ndelta36*uy3;
		Ndeltav[3] = Ndelta41*ux1 +  Ndelta42*uy1 +  Ndelta43*ux2 +  Ndelta44*uy2 +  Ndelta45*ux3 +  Ndelta46*uy3;
		Ndeltav[4] = Ndelta51*ux1 +  Ndelta52*uy1 +  Ndelta53*ux2 +  Ndelta54*uy2 +  Ndelta55*ux3 +  Ndelta56*uy3;
		Ndeltav[5] = Ndelta61*ux1 +  Ndelta62*uy1 +  Ndelta63*ux2 +  Ndelta64*uy2 +  Ndelta65*ux3 +  Ndelta66*uy3;

		Gdeltav[0] = Gdelta11*p1 + Gdelta12*p2 + Gdelta13*p3;
		Gdeltav[1] = Gdelta21*p1 + Gdelta22*p2 + Gdelta23*p3;
		Gdeltav[2] = Gdelta31*p1 + Gdelta32*p2 + Gdelta33*p3;
		Gdeltav[3] = Gdelta41*p1 + Gdelta42*p2 + Gdelta43*p3;
		Gdeltav[4] = Gdelta51*p1 + Gdelta52*p2 + Gdelta53*p3;
		Gdeltav[5] = Gdelta61*p1 + Gdelta62*p2 + Gdelta63*p3;

		Nphiv[0] = Nphi11*ux1 + Nphi12*uy1 + Nphi13*ux2 + Nphi14*uy2 + Nphi15*ux3 + Nphi16*uy3;
		Nphiv[1] = Nphi21*ux1 + Nphi22*uy1 + Nphi23*ux2 + Nphi24*uy2 + Nphi25*ux3 + Nphi26*uy3;
		Nphiv[2] = Nphi31*ux1 + Nphi32*uy1 + Nphi33*ux2 + Nphi34*uy2 + Nphi35*ux3 + Nphi36*uy3;

		Gphiv[0] = Gphi11*p1 + Gphi12*p2 + Gphi13*p3; 
		Gphiv[1] = Gphi21*p1 + Gphi22*p2 + Gphi23*p3;
		Gphiv[2] = Gphi31*p1 + Gphi32*p2 + Gphi33*p3;

		//*****************Residuo do metodo de Newton*****************
		
		Res[0] = Fe[0] - (Nv[0] - Ndeltav[0] + Kv[0] + Ksv[0] - ( Gv[0] - Gdeltav[0] ) );
		Res[1] = Fe[1] - (Nv[1] - Ndeltav[1] + Kv[1] + Ksv[1] - ( Gv[1] - Gdeltav[1] ) );
		Res[2] = Fe[2] - ( GTv[0] - Nphiv[0] + Gphiv[0] );
		Res[3] = Fe[3] - (Nv[2] - Ndeltav[2] + Kv[2] + Ksv[2] - ( Gv[2] - Gdeltav[2] ) );
		Res[4] = Fe[4] - (Nv[3] - Ndeltav[3] + Kv[3] + Ksv[3] - ( Gv[3] - Gdeltav[3] ) );
		Res[5] = Fe[5] - ( GTv[1] - Nphiv[1] + Gphiv[1] );
		Res[6] = Fe[6] - (Nv[4] - Ndeltav[4] + Kv[4] + Ksv[4] - ( Gv[4] - Gdeltav[4] ) );
		Res[7] = Fe[7] - (Nv[5] - Ndeltav[5] + Kv[5] + Ksv[5] - ( Gv[5] - Gdeltav[5] ) );
		Res[8] = Fe[8] - ( GTv[2] - Nphiv[2] + Gphiv[2] );
				
		//**********************Tangent matrix*******************************
				
		//if(varglobal < 50000){ //varglobal (<5 ISI, Ponto Fixo) (>=5 NI. Met. de Newton)

		Ke[0][0] = N11 - Ndelta11 + K11 + Ks11;
		Ke[0][1] =       K12 + Ks12;
		Ke[0][2] = -( GT11 - Gdelta11 );
		Ke[0][3] = N13 - Ndelta13 + K13 + Ks13;
		Ke[0][4] =		  K14 + Ks14;
		Ke[0][5] = -( GT11 - Gdelta12 );
		Ke[0][6] = N15 - Ndelta15 + K15 + Ks15;
		Ke[0][7] = 		  K16 + Ks16;
		Ke[0][8] = -( GT11 - Gdelta13 );
	
		Ke[1][0] = 		  K12 + Ks21;
		Ke[1][1] = N11 - Ndelta22 + K22 + Ks22;
		Ke[1][2] = -( GT12 - Gdelta21 );
		Ke[1][3] = 		  K23 + Ks23;
		Ke[1][4] = N13 - Ndelta24 + K24 + Ks24;
		Ke[1][5] = -( GT12 - Gdelta22 );
		Ke[1][6] = 		  K25 + Ks25;
		Ke[1][7] = N15 - Ndelta26 + K26 + Ks26;
		Ke[1][8] = -( GT12 - Gdelta23 );
		
		Ke[2][0] = GT11 - Gdelta11;
		Ke[2][1] = GT12 - Gdelta21;
		Ke[2][2] = Gphi11;
		Ke[2][3] = GT13 - Gdelta31;
		Ke[2][4] = GT14 - Gdelta41;
		Ke[2][5] = Gphi12;
		Ke[2][6] = GT15 - Gdelta51;
		Ke[2][7] = GT16 - Gdelta61;
		Ke[2][8] = Gphi13;

		Ke[3][0] = N11 - Ndelta13 + K13 + Ks31;
		Ke[3][1] =        K23 + Ks32;
		Ke[3][2] = -( GT13 - Gdelta31 );
		Ke[3][3] = N13 - Ndelta33 + K33 + Ks33;
		Ke[3][4] =	  K34 + Ks34;
		Ke[3][5] = -( GT13 - Gdelta32 );
		Ke[3][6] = N15 - Ndelta35 + K35 + Ks35;
		Ke[3][7] = 	  K36 + Ks36;
		Ke[3][8] = -( GT13 - Gdelta33 );

		Ke[4][0] = 	 K14 + Ks41;
		Ke[4][1] = N11 - Ndelta24 + K24 + Ks42;
		Ke[4][2] = -( GT14 - Gdelta41 );
		Ke[4][3] = 	 K34 + Ks43;
		Ke[4][4] = N13 - Ndelta44 + K44 + Ks44;
		Ke[4][5] = -( GT14 - Gdelta42 );
		Ke[4][6] = 	 K45 + Ks45;
		Ke[4][7] = N15 - Ndelta46 + K46 + Ks46;
		Ke[4][8] = -( GT14 - Gdelta43 );

		Ke[5][0] = GT11 - Gdelta12;
		Ke[5][1] = GT12 - Gdelta22;
		Ke[5][2] = Gphi12;
		Ke[5][3] = GT13 - Gdelta32;
		Ke[5][4] = GT14 - Gdelta42;
		Ke[5][5] = Gphi22;
		Ke[5][6] = GT15 - Gdelta52;
		Ke[5][7] = GT16 - Gdelta62;
		Ke[5][8] = Gphi23;

		Ke[6][0] = N11 - Ndelta15 + K15 + Ks51;
		Ke[6][1] =       K25 + Ks52;
		Ke[6][2] = -( GT15 - Gdelta51 );
		Ke[6][3] = N13 + Ndelta35 + K35 + Ks53;
		Ke[6][4] =	 K45 + Ks54;
		Ke[6][5] = -( GT15 - Gdelta52 );
		Ke[6][6] = N15 - Ndelta55 + K55 + Ks55;
		Ke[6][7] = 	 K56 + Ks56;
		Ke[6][8] = -( GT15 - Gdelta53 );

		Ke[7][0] = 	 K16 + Ks61;
		Ke[7][1] = N11 - Ndelta26 + K26 + Ks62;
		Ke[7][2] = -( GT16 - Gdelta61 );
		Ke[7][3] = 	 K36 + Ks63;
		Ke[7][4] = N13 - Ndelta46 + K46 + Ks64;
		Ke[7][5] = -( GT16 - Gdelta62 );
		Ke[7][6] = 	 K56 + Ks65;
		Ke[7][7] = N15 - Ndelta66 + K66 + Ks66;
		Ke[7][8] = -( GT16 - Gdelta63 );

		Ke[8][0] = GT11 - Gdelta13;
		Ke[8][1] = GT12 - Gdelta23;
		Ke[8][2] = Gphi13;
		Ke[8][3] = GT13 - Gdelta33;
		Ke[8][4] = GT14 - Gdelta43;
		Ke[8][5] = Gphi23;
		Ke[8][6] = GT15 - Gdelta53;
		Ke[8][7] = GT16 - Gdelta63;
		Ke[8][8] = Gphi33;
		
		//}



		/* direciona para um lixo durante o produto matriz vetor os elementos da matriz
		 que deveriam ser nulos (os que foram para o vetor força devido a c.c.). 
		Nos da a posicao de assemble global do vetor forca
		Aqui neq == numequacao + 1. */
/*		for (i = 0; i<NDOF; i++){
			if (Node[J1].id[i]==-1)
				lm_aux[i] = neq;
			else
				lm_aux[i] = Node[J1].id[i];

			if (Node[J2].id[i]==-1)
				lm_aux[i + NDOF] = neq;
			else
				lm_aux[i + NDOF] = Node[J2].id[i];
		
			if (Node[J3].id[i]==-1)
				lm_aux[i + 2*NDOF] = neq;
			else
				lm_aux[i + 2*NDOF] = Node[J3].id[i];
		}
*/
		// Assemble global do vetor independente F de Au=F 
		for (i = 0; i < 9; i++)
			F[lm[E][i]] +=  Res[i];
		
		F[NEQ] = 0;


		
		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, E, Ke);
		
	}//for elemento

	free(U);		
	varglobal++;

	return 0;

}
