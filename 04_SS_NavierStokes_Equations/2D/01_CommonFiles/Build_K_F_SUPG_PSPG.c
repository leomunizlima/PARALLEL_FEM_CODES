#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

static int varglobal = 0;

int Build_K_F_SUPG_PSPG(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int i, NEQ, nel;//, op;
	int E, J1, J2, J3;
	//double Uref, Lref;
	double Re, x, y;
	double TwoA, invArea, Area, tau, h, visc, rho, invrho, nu, unorm, auxtau, tau_l; 
	double third=1.0/3.0, sixth = 1.0/6.0, twelfth = 1.0/12.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Ke[9][9], Fe[9], ux1, ux2, ux3, uy1, uy2, uy3, ux, uy, duxdx, duxdy, duydx, duydy;
	double fx1, fx2, fx3, fy1, fy2, fy3, fxB, fyB, f1, f2, fdelta1, fdelta2, fdelta3, fdelta4, fdelta5, fdelta6, fphi1, fphi2, fphi3;
	double C1, C2, C3;
	double Nv[6], Ndeltav[6], Kv[6], Ksv[6], Gv[6], Gdeltav[6], GTv[3], Nphiv[3], Gphiv[3], Res[9];
	//double M11, M13, Mdelta11, Mdelta31, Mdelta51;
	double K11, K12, K13, K14, K15, K16, K22, K23, K24, K25, K26, K33, K34, K35, K36, K44, K45, K46, K55, K56, K66; 
	double Ks11, Ks12, Ks13, Ks14, Ks15, Ks16, Ks21, Ks22, Ks23, Ks24, Ks25, Ks26, Ks31, Ks32, Ks33, Ks34, Ks35, Ks36;
	double Ks41, Ks42, Ks43, Ks44, Ks45, Ks46, Ks51, Ks52, Ks53, Ks54, Ks55, Ks56, Ks61, Ks62, Ks63, Ks64, Ks65, Ks66; 
	double GT11, GT12, GT13, GT14, GT15, GT16, N11, N13, N15;
	double Ndelta11, Ndelta13, Ndelta15, Ndelta22, Ndelta24, Ndelta26, Ndelta33, Ndelta35, Ndelta44, Ndelta46, Ndelta55, Ndelta66;
	double Gdelta11, Gdelta12, Gdelta13, Gdelta21, Gdelta22, Gdelta23, Gdelta31, Gdelta32, Gdelta33, Gdelta41, Gdelta42, Gdelta43;
	double Gdelta51, Gdelta52, Gdelta53, Gdelta61, Gdelta62, Gdelta63;
	//double Mphi11, Mphi12, Mphi21, Mphi22, Mphi31, Mphi32;
	double Gphi11, Gphi12, Gphi13, Gphi22, Gphi23, Gphi33;
	double Np11, Np12, Np21, Np22, Np31, Np32, Np41, Np42;
	double Npdelta11, Npdelta12, Npdelta21, Npdelta22, Npdelta31, Npdelta32, Npdelta41, Npdelta42, Npdelta51, Npdelta52, Npdelta61, Npdelta62; 
	double Nppdelta11, Nppdelta12, Nppdelta21, Nppdelta22, Nppdelta31, Nppdelta32, Nppdelta41, Nppdelta42, Nppdelta51, Nppdelta52, Nppdelta61, Nppdelta62;
	double Npphi11, Npphi12, Npphi21, Npphi22, Npphi31, Npphi32;
	double p1, p2, p3, p, dpdx, dpdy, duxduy;
	
	double *F = FemStructs->F;
	int **lm = FemStructs->lm;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	nel = Parameters->nel;
	NEQ = Parameters->NEQ;
	
	dzero(NEQ+1, F);
	setzeros(Parameters,MatrixData);

	double *U = (double*) mycalloc("U of 'Build_K_F_SUPG_PSPG'", 3*Parameters->nnodes, sizeof(double));
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
		ux = third*( ux1 + ux2 + ux3 );
		uy = third*( uy1 + uy2 + uy3 );
		
		//*****************************Gradient of u********************
		duxdx = 0.5*invArea*( ux1*y23 + ux2*y31 + ux3*y12 ); 
		duxdy = 0.5*invArea*( ux1*x32 + ux2*x13 + ux3*x21 );
		duydx = 0.5*invArea*( uy1*y23 + uy2*y31 + uy3*y12 );
		duydy = 0.5*invArea*( uy1*x32 + uy2*x13 + uy3*x21 );
				
		duxduy = ux1*y23 + uy1*x32 + ux2*y31 + uy2*x13 + ux3*y12 + uy3*x21;

		//*****************p at barycentre of the element***************
		p = third*( p1 + p2 + p3 );

		//**************************Gradient of p**********************
		dpdx = 0.5*invArea*( p1*y23 + p2*y31 + p3*y12 );
		dpdy = 0.5*invArea*( p1*x32 + p2*x13 + p3*x21 );

		//*****************Coefficients C_i******************************
		C1 = 0.5*invArea*( ux*y23 + uy*x32 );
		C2 = 0.5*invArea*( ux*y31 + uy*x13 );
		C3 = 0.5*invArea*( ux*y12 + uy*x21 );

		//*************Calculation of tau_SUPG=tau_PSPG=tau*********
		Re = Parameters->ReynoldsNumber;		
		h = sqrt( 4*Area/PI );		
		visc =1./Re;
		rho = 1.;		
		//visc = 1./sqrt(Re);
		//rho = sqrt(Re);		

		//visc = (*Viscosity)(); 
		//rho = (*Rho)();
		//Uref = (*U_ref)();
		//Lref = (*L_ref)();
		invrho = 1.0/rho;
		nu = visc/rho;
		//h = sqrt( 4*Area/PI );		
		unorm = sqrt( ux*ux + uy*uy );
		auxtau = ( 2*unorm/h )*( 2*unorm/h ) + 9*( 4*nu/(h*h) )*( 4*nu/(h*h) );
		auxtau = sqrt( auxtau);
		tau = 1.0/auxtau;
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
		
		//***************Matrices of the SUPG formulation****************

		//*******************Matrix Mdelta=tau*N^T***********************
		//Mdelta11 = tau*N11;
		//Mdelta31 = tau*N13;
		//Mdelta51 = tau*N15;

		//**********************Matrix Ndelta****************************
		Ndelta11 = rho*tau*Area*C1*C1;
		Ndelta13 = rho*tau*Area*C1*C2;
		Ndelta15 = rho*tau*Area*C1*C3;
	
		Ndelta22 = Ndelta11;
		Ndelta24 = Ndelta13;
		Ndelta26 = Ndelta15;
	
		Ndelta33 = rho*tau*Area*C2*C2;
		Ndelta35 = rho*tau*Area*C2*C3;
	
		Ndelta44 = Ndelta33;
		Ndelta46 = Ndelta35;
	
		Ndelta55 = rho*tau*Area*C3*C3;
	
		Ndelta66 = Ndelta55;
		
		//**********************Matrix Gdelta****************************
		Gdelta11 = 0.5*tau*C1*y23;
		Gdelta12 = 0.5*tau*C1*y31;
		Gdelta13 = 0.5*tau*C1*y12;

		Gdelta21 = 0.5*tau*C1*x32;
		Gdelta22 = 0.5*tau*C1*x13;
		Gdelta23 = 0.5*tau*C1*x21;

		Gdelta31 = 0.5*tau*C2*y23;
		Gdelta32 = 0.5*tau*C2*y31;
		Gdelta33 = 0.5*tau*C2*y12;

		Gdelta41 = 0.5*tau*C2*x32;
		Gdelta42 = 0.5*tau*C2*x13;
		Gdelta43 = 0.5*tau*C2*x21;

		Gdelta51 = 0.5*tau*C3*y23;
		Gdelta52 = 0.5*tau*C3*y31;
		Gdelta53 = 0.5*tau*C3*y12;

		Gdelta61 = 0.5*tau*C3*x32;
		Gdelta62 = 0.5*tau*C3*x13;
		Gdelta63 = 0.5*tau*C3*x21;

		//***************Matrices of the PSPG formulation****************

		//**********************Matrix Mphi******************************
		//Mphi11 = tau*sixth*y23;
		//Mphi12 = tau*sixth*x32;
	
		//Mphi21 = tau*sixth*y31;
		//Mphi22 = tau*sixth*x13;
	
		//Mphi31 = tau*sixth*y12;
		//Mphi32 = tau*sixth*x21;

		//***************** Matrix Nphi = Gdelta^T ***********************

		//**********************Matrix Gphi*******************************
		Gphi11 = invrho*0.25*tau*invArea*( y23*y23 + x32*x32 );
		Gphi12 = invrho*0.25*tau*invArea*( y23*y31 + x32*x13 );
		Gphi13 = invrho*0.25*tau*invArea*( y23*y12 + x32*x21 );
		
		Gphi22 = invrho*0.25*tau*invArea*( y31*y31 + x13*x13 );
		Gphi23 = invrho*0.25*tau*invArea*( y31*y12 + x13*x21 );
		
		Gphi33 = invrho*0.25*tau*invArea*( y12*y12 + x21*x21 );	
		
		//*********************Incremental matrices***********************
		
		//*************************Matrix Np*******************************
		Np11 = rho*Area*sixth*duxdx;
		Np12 = rho*Area*sixth*duxdy;
		Np21 = rho*Area*sixth*duydx;
		Np22 = rho*Area*sixth*duydy;

		Np31 = rho*Area*twelfth*duxdx;
		Np32 = rho*Area*twelfth*duxdy;
		Np41 = rho*Area*twelfth*duydx;
		Np42 = rho*Area*twelfth*duydy;
		
		//**********************Matrix Npdelta*******************************
		Npdelta11 = rho*tau*Area*third*C1*duxdx;
		Npdelta12 = rho*tau*Area*third*C1*duxdy;
		Npdelta21 = rho*tau*Area*third*C1*duydx;
		Npdelta22 = rho*tau*Area*third*C1*duydy;
	
		Npdelta31 = rho*tau*Area*third*C2*duxdx;
		Npdelta32 = rho*tau*Area*third*C2*duxdy;
		Npdelta41 = rho*tau*Area*third*C2*duydx;
		Npdelta42 = rho*tau*Area*third*C2*duydy;
	
		Npdelta51 = rho*tau*Area*third*C3*duxdx;
		Npdelta52 = rho*tau*Area*third*C3*duxdy;
		Npdelta61 = rho*tau*Area*third*C3*duydx;
		Npdelta62 = rho*tau*Area*third*C3*duydy;

		//*********************Matrix Nppdelta*******************************
		Nppdelta11 = rho*tau*sixth*y23*( duxdx*ux + duxdy*uy );
		Nppdelta12 = rho*tau*sixth*x32*( duxdx*ux + duxdy*uy );
		Nppdelta21 = rho*tau*sixth*y23*( duydx*ux + duydy*uy );
		Nppdelta22 = rho*tau*sixth*x32*( duydx*ux + duydy*uy );

		Nppdelta31 = rho*tau*sixth*y31*( duxdx*ux + duxdy*uy );
		Nppdelta32 = rho*tau*sixth*x13*( duxdx*ux + duxdy*uy );
		Nppdelta41 = rho*tau*sixth*y31*( duydx*ux + duydy*uy );
		Nppdelta42 = rho*tau*sixth*x13*( duydx*ux + duydy*uy );

		Nppdelta51 = rho*tau*sixth*y12*( duxdx*ux + duxdy*uy );
		Nppdelta52 = rho*tau*sixth*x21*( duxdx*ux + duxdy*uy );
		Nppdelta61 = rho*tau*sixth*y12*( duydx*ux + duydy*uy );
		Nppdelta62 = rho*tau*sixth*x21*( duydx*ux + duydy*uy );

		//************************Matrix Npphi*******************************
		Npphi11 = tau*sixth*( duxdx*y23 + duydx*x32 );
		Npphi12 = tau*sixth*( duxdy*y23 + duydy*x32 );

		Npphi21 = tau*sixth*( duxdx*y31 + duydx*x13 );
		Npphi22 = tau*sixth*( duxdy*y31 + duydy*x13 );

		Npphi31 = tau*sixth*( duxdx*y12 + duydx*x21 );
		Npphi32 = tau*sixth*( duxdy*y12 + duydy*x21 );
	
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


		Ndeltav[0] = rho*tau*Area*C1*( duxdx*ux + duxdy*uy );
		Ndeltav[1] = rho*tau*Area*C1*( duydx*ux + duydy*uy );
		Ndeltav[2] = rho*tau*Area*C2*( duxdx*ux + duxdy*uy );
		Ndeltav[3] = rho*tau*Area*C2*( duydx*ux + duydy*uy );
		Ndeltav[4] = rho*tau*Area*C3*( duxdx*ux + duxdy*uy );
		Ndeltav[5] = rho*tau*Area*C3*( duydx*ux + duydy*uy );


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

		Gdeltav[0] = tau*Area*C1*dpdx;
		Gdeltav[1] = tau*Area*C1*dpdy;
		Gdeltav[2] = tau*Area*C2*dpdx;
		Gdeltav[3] = tau*Area*C2*dpdy;
		Gdeltav[4] = tau*Area*C3*dpdx;
		Gdeltav[5] = tau*Area*C3*dpdy;
		
		GTv[0] = third*Area*( duxdx + duydy);
		GTv[1] = third*Area*( duxdx + duydy);
		GTv[2] = third*Area*( duxdx + duydy);
		
		Nphiv[0] = 0.5*tau*(( duxdx*ux + duxdy*uy )*y23 + ( duydx*ux + duydy*uy )*x32);
		Nphiv[1] = 0.5*tau*(( duxdx*ux + duxdy*uy )*y31 + ( duydx*ux + duydy*uy )*x13);
		Nphiv[2] = 0.5*tau*(( duxdx*ux + duxdy*uy )*y12 + ( duydx*ux + duydy*uy )*x21);

		Gphiv[0] = 0.5*tau*invrho*( y23*dpdx + x32*dpdy); 
		Gphiv[1] = 0.5*tau*invrho*( y31*dpdx + x13*dpdy);
		Gphiv[2] = 0.5*tau*invrho*( y12*dpdx + x21*dpdy);

		//Residue 
		
		Res[0] = Fe[0] - (Nv[0] + Ndeltav[0] + Kv[0] + Ksv[0] - ( Gv[0] + Gdeltav[0] ) );
		Res[1] = Fe[1] - (Nv[1] + Ndeltav[1] + Kv[1] + Ksv[1] - ( Gv[1] + Gdeltav[1] ) );
		Res[2] = Fe[2] - ( GTv[0] + Nphiv[0] + Gphiv[0] );
		Res[3] = Fe[3] - (Nv[2] + Ndeltav[2] + Kv[2] + Ksv[2] - ( Gv[2] + Gdeltav[2] ) );
		Res[4] = Fe[4] - (Nv[3] + Ndeltav[3] + Kv[3] + Ksv[3] - ( Gv[3] + Gdeltav[3] ) );
		Res[5] = Fe[5] - ( GTv[1] + Nphiv[1] + Gphiv[1] );
		Res[6] = Fe[6] - (Nv[4] + Ndeltav[4] + Kv[4] + Ksv[4] - ( Gv[4] + Gdeltav[4] ) );
		Res[7] = Fe[7] - (Nv[5] + Ndeltav[5] + Kv[5] + Ksv[5] - ( Gv[5] + Gdeltav[5] ) );
		Res[8] = Fe[8] - ( GTv[2] + Nphiv[2] + Gphiv[2] );
				
		//**********************Tangent matrix*******************************
				
		if(varglobal < 5){ //varglobal (<5 ISI, Ponto Fixo) (>=5 NI. Met. de Newton)

			Ke[0][0] = N11 + Ndelta11 + K11 + Ks11;
			Ke[0][1] =       K12 + Ks12;
			Ke[0][2] = -( GT11 + Gdelta11 );
			Ke[0][3] = N13 + Ndelta13 + K13 + Ks13;
			Ke[0][4] =		  K14 + Ks14;
			Ke[0][5] = -( GT11 + Gdelta12 );
			Ke[0][6] = N15 + Ndelta15 + K15 + Ks15;
			Ke[0][7] = 		  K16 + Ks16;
			Ke[0][8] = -( GT11 + Gdelta13 );
	
			Ke[1][0] = 		  K12 + Ks21;
			Ke[1][1] = N11 + Ndelta22 + K22 + Ks22;
			Ke[1][2] = -( GT12 + Gdelta21 );
			Ke[1][3] = 		  K23 + Ks23;
			Ke[1][4] = N13 + Ndelta24 + K24 + Ks24;
			Ke[1][5] = -( GT12 + Gdelta22 );
			Ke[1][6] = 		  K25 + Ks25;
			Ke[1][7] = N15 + Ndelta26 + K26 + Ks26;
			Ke[1][8] = -( GT12 + Gdelta23 );
		
			Ke[2][0] = GT11 + Gdelta11;
			Ke[2][1] = GT12 + Gdelta21;
			Ke[2][2] = Gphi11;
			Ke[2][3] = GT13 + Gdelta31;
			Ke[2][4] = GT14 + Gdelta41;
			Ke[2][5] = Gphi12;
			Ke[2][6] = GT15 + Gdelta51;
			Ke[2][7] = GT16 + Gdelta61;
			Ke[2][8] = Gphi13;

			Ke[3][0] = N11 + Ndelta13 + K13 + Ks31;
			Ke[3][1] =        K23 + Ks32;
			Ke[3][2] = -( GT13 + Gdelta31 );
			Ke[3][3] = N13 + Ndelta33 + K33 + Ks33;
			Ke[3][4] =		   K34 + Ks34;
			Ke[3][5] = -( GT13 + Gdelta32 );
			Ke[3][6] = N15 + Ndelta35 + K35 + Ks35;
			Ke[3][7] = 		   K36 + Ks36;
			Ke[3][8] = -( GT13 + Gdelta33 );
	
			Ke[4][0] = 		  K14 + Ks41;
			Ke[4][1] = N11 + Ndelta24 + K24 + Ks42;
			Ke[4][2] = -( GT14 + Gdelta41 );
			Ke[4][3] = 		  K34 + Ks43;
			Ke[4][4] = N13 + Ndelta44 + K44 + Ks44;
			Ke[4][5] = -( GT14 + Gdelta42 );
			Ke[4][6] = 		  K45 + Ks45;
			Ke[4][7] = N15 + Ndelta46 + K46 + Ks46;
			Ke[4][8] = -( GT14 + Gdelta43 );
	
			Ke[5][0] = GT11 + Gdelta12;
			Ke[5][1] = GT12 + Gdelta22;
			Ke[5][2] = Gphi12;
			Ke[5][3] = GT13 + Gdelta32;
			Ke[5][4] = GT14 + Gdelta42;
			Ke[5][5] = Gphi22;
			Ke[5][6] = GT15 + Gdelta52;
			Ke[5][7] = GT16 + Gdelta62;
			Ke[5][8] = Gphi23;
	
			Ke[6][0] = N11 + Ndelta15 + K15 + Ks51;
			Ke[6][1] =       K25 + Ks52;
			Ke[6][2] = -( GT15 + Gdelta51 );
			Ke[6][3] = N13 + Ndelta35 + K35 + Ks53;
			Ke[6][4] =		 K45 + Ks54;
			Ke[6][5] = -( GT15 + Gdelta52 );
			Ke[6][6] = N15 + Ndelta55 + K55 + Ks55;
			Ke[6][7] = 		  K56 + Ks56;
			Ke[6][8] = -( GT15 + Gdelta53 );

			Ke[7][0] = 		  K16 + Ks61;
			Ke[7][1] = N11 + Ndelta26 + K26 + Ks62;
			Ke[7][2] = -( GT16 + Gdelta61 );
			Ke[7][3] = 		  K36 + Ks63;
			Ke[7][4] = N13 + Ndelta46 + K46 + Ks64;
			Ke[7][5] = -( GT16 + Gdelta62 );
			Ke[7][6] = 		  K56 + Ks65;
			Ke[7][7] = N15 + Ndelta66 + K66 + Ks66;
			Ke[7][8] = -( GT16 + Gdelta63 );

			Ke[8][0] = GT11 + Gdelta13;
			Ke[8][1] = GT12 + Gdelta23;
			Ke[8][2] = Gphi13;
			Ke[8][3] = GT13 + Gdelta33;
			Ke[8][4] = GT14 + Gdelta43;
			Ke[8][5] = Gphi23;
			Ke[8][6] = GT15 + Gdelta53;
			Ke[8][7] = GT16 + Gdelta63;
			Ke[8][8] = Gphi33;
		
		}else{

			Ke[0][0] = N11 + Np11 + Ndelta11 + Npdelta11 + Nppdelta11 + K11  + Ks11;
			Ke[0][1] =       Np12 			 + Npdelta12 + Nppdelta12 + K12 + Ks12;
			Ke[0][2] = -( GT11 + Gdelta11 );
			Ke[0][3] = N13 + Np31 + Ndelta13 + Npdelta11 + Nppdelta11 + K13 + Ks13;
			Ke[0][4] =		 Np32 	 		 + Npdelta12 + Nppdelta12 + K14 + Ks14;
			Ke[0][5] = -( GT11 + Gdelta12 );
			Ke[0][6] = N15 + Np31 + Ndelta15 + Npdelta11 + Nppdelta11 + K15 + Ks15;
			Ke[0][7] = 		 Np32 	 		 + Npdelta12 + Nppdelta12 + K16 + Ks16;
			Ke[0][8] = -( GT11 + Gdelta13 );
	
			Ke[1][0] = 		 Np21            + Npdelta21 + Nppdelta21 + K12 + Ks21;
			Ke[1][1] = N11 + Np22 + Ndelta22 + Npdelta22 + Nppdelta22 + K22 + Ks22;
			Ke[1][2] = -( GT12 + Gdelta21 );
			Ke[1][3] = 		 Np41			 + Npdelta21 + Nppdelta21 + K23 + Ks23;
			Ke[1][4] = N13 + Np42 + Ndelta24 + Npdelta22 + Nppdelta22 + K24 + Ks24;
			Ke[1][5] = -( GT12 + Gdelta22 );
			Ke[1][6] = 		 Np41 			 + Npdelta21 + Nppdelta21 + K25 + Ks25;
			Ke[1][7] = N15 + Np42 + Ndelta26 + Npdelta22 + Nppdelta22 + K26 + Ks26;
			Ke[1][8] = -( GT12 + Gdelta23 );
		
			Ke[2][0] = GT11 + Gdelta11 + Npphi11;
			Ke[2][1] = GT12 + Gdelta21 + Npphi12;
			Ke[2][2] = Gphi11;
			Ke[2][3] = GT13 + Gdelta31 + Npphi11;
			Ke[2][4] = GT14 + Gdelta41 + Npphi12;
			Ke[2][5] = Gphi12;
			Ke[2][6] = GT15 + Gdelta51 + Npphi11;
			Ke[2][7] = GT16 + Gdelta61 + Npphi12;
			Ke[2][8] = Gphi13;
	
			Ke[3][0] = N11 + Np31 + Ndelta13 + Npdelta31 + Nppdelta31 + K13 + Ks31;
			Ke[3][1] =       Np32 			 + Npdelta32 + Nppdelta32 + K23 + Ks32;
			Ke[3][2] = -( GT13 + Gdelta31 );
			Ke[3][3] = N13 + Np11 + Ndelta33 + Npdelta31 + Nppdelta31 + K33 + Ks33;
			Ke[3][4] =		 Np12 	 		 + Npdelta32 + Nppdelta32 + K34 + Ks34;
			Ke[3][5] = -( GT13 + Gdelta32 );
			Ke[3][6] = N15 + Np31 + Ndelta35 + Npdelta31 + Nppdelta31 + K35 + Ks35;
			Ke[3][7] = 		 Np32 	 		 + Npdelta32 + Nppdelta32 + K36 + Ks36;
			Ke[3][8] = -( GT13 + Gdelta33 );
	

			Ke[4][0] = 		 Np41            + Npdelta41 + Nppdelta41 + K14 + Ks41;
			Ke[4][1] = N11 + Np42 + Ndelta24 + Npdelta42 + Nppdelta42 + K24 + Ks42;
			Ke[4][2] = -( GT14 + Gdelta41 );
			Ke[4][3] = 		 Np21			 + Npdelta41 + Nppdelta41 + K34 + Ks43;
			Ke[4][4] = N13 + Np22 + Ndelta44 + Npdelta42 + Nppdelta42 + K44 + Ks44;
			Ke[4][5] = -( GT14 + Gdelta42 );
			Ke[4][6] = 		 Np41 			 + Npdelta41 + Nppdelta41 + K45 + Ks45;
			Ke[4][7] = N15 + Np42 + Ndelta46 + Npdelta42 + Nppdelta42 + K46 + Ks46;
			Ke[4][8] = -( GT14 + Gdelta43 );
	
			Ke[5][0] = GT11 + Gdelta12 + Npphi21;
			Ke[5][1] = GT12 + Gdelta22 + Npphi22;
			Ke[5][2] = Gphi12;
			Ke[5][3] = GT13 + Gdelta32 + Npphi21;
			Ke[5][4] = GT14 + Gdelta42 + Npphi22;
			Ke[5][5] = Gphi22;
			Ke[5][6] = GT15 + Gdelta52 + Npphi21;
			Ke[5][7] = GT16 + Gdelta62 + Npphi22;
			Ke[5][8] = Gphi23;
	
			Ke[6][0] = N11 + Np31 + Ndelta15 + Npdelta51 + Nppdelta51 + K15 + Ks51;
			Ke[6][1] =       Np32 			 + Npdelta52 + Nppdelta52 + K25 + Ks52;
			Ke[6][2] = -( GT15 + Gdelta51 );
			Ke[6][3] = N13 + Np31 + Ndelta35 + Npdelta51 + Nppdelta51 + K35 + Ks53;
			Ke[6][4] =		 Np32 	 		 + Npdelta52 + Nppdelta52 + K45 + Ks54;
			Ke[6][5] = -( GT15 + Gdelta52 );
			Ke[6][6] = N15 + Np11 + Ndelta55 + Npdelta51 + Nppdelta51 + K55 + Ks55;
			Ke[6][7] = 		 Np12 	 		 + Npdelta52 + Nppdelta52 + K56 + Ks56;
			Ke[6][8] = -( GT15 + Gdelta53 );

			Ke[7][0] = 		 Np41            + Npdelta61 + Nppdelta61 + K16 + Ks61;
			Ke[7][1] = N11 + Np42 + Ndelta26 + Npdelta62 + Nppdelta62 + K26 + Ks62;
			Ke[7][2] = -( GT16 + Gdelta61 );
			Ke[7][3] = 		 Np41			 + Npdelta61 + Nppdelta61 + K36 + Ks63;
			Ke[7][4] = N13 + Np42 + Ndelta46 + Npdelta62 + Nppdelta62 + K46 + Ks64;
			Ke[7][5] = -( GT16 + Gdelta62 );
			Ke[7][6] = 		 Np21 			 + Npdelta61 + Nppdelta61 + K56 + Ks65;
			Ke[7][7] = N15 + Np22 + Ndelta66 + Npdelta62 + Nppdelta62 + K66 + Ks66;
			Ke[7][8] = -( GT16 + Gdelta63 );

			Ke[8][0] = GT11 + Gdelta13 + Npphi31;
			Ke[8][1] = GT12 + Gdelta23 + Npphi32;
			Ke[8][2] = Gphi13;
			Ke[8][3] = GT13 + Gdelta33 + Npphi31;
			Ke[8][4] = GT14 + Gdelta43 + Npphi32;
			Ke[8][5] = Gphi23;
			Ke[8][6] = GT15 + Gdelta53 + Npphi31;
			Ke[8][7] = GT16 + Gdelta63 + Npphi32;
			Ke[8][8] = Gphi33;
		}


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
