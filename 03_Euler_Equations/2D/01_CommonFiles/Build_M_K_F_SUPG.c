#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_M_K_F_SUPG(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int e, i, j;
	int J1, J2, J3, nel, NEQ;
	double x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21;
	double Area, twoArea, third = 1.0/3.0, sixth = 1.0/6.0;
	double c, h, abs_vbeta, betax, betay, betaxy, CFL;
	double delta, psi;
	double tau, tau_a, tau_d, tau_t;
	double *U, *dU;
	double gamma;
	double Me[12][12], Ke[12][12], Ub[4], dUb[4], Ue[12], dUe[12], gradUx[4], gradUy[4];
	double tolerance;
	double alpha = Parameters->Alpha_Build;
	double delta_t = Parameters->DeltaT_Build;
	double *R = FemStructs->F;
	double *delta_old = FemStructs->delta_old;
	int **lm = FemStructs->lm;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	
	tolerance = FemStructs->AuxBuild->tolerance;
	nel = Parameters->nel;
	NEQ = Parameters->NEQ;
	
	dzero(NEQ+1, R);
	setzeros(Parameters,MatrixData);

	U = (double*) mycalloc("U of 'Build_M_F_DD_Transiente'", 4*Parameters->nnodes, sizeof(double));
	dU = (double*) mycalloc("dU of 'Build_M_F_DD_Transiente'", 4*Parameters->nnodes, sizeof(double));
	eval_U_dU(Parameters,FemStructs,FemFunctions,U,dU);


	for (e = 0; e < nel; e++){
		// *** Nos do elemento
		J1 = Element[e].Vertex[0];
		J2 = Element[e].Vertex[1];
		J3 = Element[e].Vertex[2];

		// *** Coordenadas nodais e operador diferencial
		x1 = Node[J1].x;
		x2 = Node[J2].x;
		x3 = Node[J3].x;
		
		y1 = Node[J1].y;
		y2 = Node[J2].y;
		y3 = Node[J3].y;
	
		gamma = FemFunctions->gamma(third*(x1+x2+x3), third*(y1+y2+y3));

		y23 = y2 - y3;
		y31 = y3 - y1;
		y12 = y1 - y2;

		x32 = x3 - x2;
		x13 = x1 - x3;
		x21 = x2 - x1;

		// *** Area do elemento
		twoArea = fabs(x21 * y31 - x13 * y12);
		Area = 0.5*twoArea; 

		// Fill ' Ue ' with the values ​​of the previous solution or workaround for each degree of liberty of the element node
		//Density
		Ue[0] = U[4*J1];
		Ue[4] = U[4*J2];
		Ue[8] = U[4*J3];

		dUe[0] = dU[4*J1];
		dUe[4] = dU[4*J2];
		dUe[8] = dU[4*J3];

		//x Velocity
		Ue[1] = U[4*J1+1];
		Ue[5] = U[4*J2+1];
		Ue[9] = U[4*J3+1];

		dUe[1] = dU[4*J1+1];
		dUe[5] = dU[4*J2+1];
		dUe[9] = dU[4*J3+1];

		//y Velocity
		Ue[2] = U[4*J1+2];
		Ue[6] = U[4*J2+2];
		Ue[10] = U[4*J3+2];

		dUe[2] = dU[4*J1+2];
		dUe[6] = dU[4*J2+2];
		dUe[10] = dU[4*J3+2];

		//Energy
		Ue[3] = U[4*J1+3];
		Ue[7] = U[4*J2+3];
		Ue[11] = U[4*J3+3];

		dUe[3] = dU[4*J1+3];
		dUe[7] = dU[4*J2+3];
		dUe[11] = dU[4*J3+3];

		// Baricentro do triangulo
		for (i = 0; i < 4; i++){
			Ub[i] = (Ue[i] + Ue[i+4] + Ue[i+8]) * third;
			dUb[i] = (dUe[i] + dUe[i+4] + dUe[i+8]) * third;
		}

		// *** Calculo do Gradiente (gradu = Bu)
		for (i = 0; i < 4; i++)
		{
			gradUx[i] = (Ue[i]*y23 + Ue[i+4]*y31 + Ue[i+8]*y12) / twoArea;
			gradUy[i] = (Ue[i]*x32 + Ue[i+4]*x13 + Ue[i+8]*x21) / twoArea;
		}

		// *** Matrizes Jacobianas Ax e Ay: variaveis primitivas
		double v1 = Ub[1] / Ub[0];
		double v2 = Ub[2] / Ub[0];
		double E  = Ub[3] / Ub[0];
		double norma_U23 = (Ub[1] * Ub[1]) + (Ub[2] * Ub[2]);

		// *** Jacobian matrices Ax and Ay calculations
		double Ax[4][4];
		double Ay[4][4];
	
		FemFunctions->Ax_Ay_calculations(gamma,Parameters->Mach,Ub, Ax, Ay);

		// *** Inversa da Metrica Riemmaniana A0^(-1)
		double rho_i = Ub[3] - ((norma_U23 * 0.5) / Ub[0]); 
		
		// Variaveis de entropia do baricentro
		double V2, V3, V4;
		V2 = Ub[1] / rho_i;
		V3 = Ub[2] / rho_i;
		V4 = (- Ub[0]) / rho_i;

		// Coeficientes de A0^(-1)
		double rhoi_e    = - 1.0 / (rho_i * V4);
		double k1        = (V2*V2 + V3*V3) * 0.5 / V4;
		double A0[4][4];

		A0[0][0] = rhoi_e * (k1 * k1 + gamma);
		A0[0][1] = rhoi_e * (k1 * V2);
		A0[0][2] = rhoi_e * (k1 * V3);
		A0[0][3] = rhoi_e * ((k1 + 1.0) * V4);
		A0[1][0] = A0[0][1];
		A0[1][1] = rhoi_e * (V2 * V2 - V4);
		A0[1][2] = rhoi_e * (V2 * V3);
		A0[1][3] = rhoi_e * (V2 * V4);
		A0[2][0] = A0[0][2];
		A0[2][1] = A0[1][2];
		A0[2][2] = rhoi_e * (V3 * V3 - V4);
		A0[2][3] = rhoi_e * (V3 * V4);
		A0[3][0] = A0[0][3];
		A0[3][1] = A0[1][3];
		A0[3][2] = A0[2][3];
		A0[3][3] = rhoi_e * V4 * V4;      

		// *** Coeficientes dos produtos AiAj, onde i,j=x,y

		double Axx11, Axx12, Axx13, Axx14, Axx21, Axx22, Axx23, Axx24, Axx31, Axx32, Axx33, Axx34, Axx41, Axx42, Axx43, Axx44;
		// calcula Axx = AxAx
		Axx11 = Ax[1][0];
		Axx12 = Ax[1][1];
		Axx13 = Ax[1][2];
		Axx14 = Ax[1][3];
		Axx21 = Ax[1][1]*Ax[1][0] + Ax[1][2]*Ax[2][0] + Ax[1][3]*Ax[3][0];
		Axx22 =      Ax[1][0]     + Ax[1][1]*Ax[1][1] + Ax[1][2]*Ax[2][1] + Ax[1][3]*Ax[3][1];
		Axx23 = Ax[1][1]*Ax[1][2] + Ax[1][2]*Ax[2][2] + Ax[1][3]*Ax[3][2];
		Axx24 = Ax[1][1]*Ax[1][3] + Ax[1][3]*Ax[3][3];
		Axx31 = Ax[2][1]*Ax[1][0] + Ax[2][2]*Ax[2][0];
		Axx32 =      Ax[2][0]     + Ax[2][1]*Ax[1][1] + Ax[2][2]*Ax[2][1];
		Axx33 = Ax[2][1]*Ax[1][2] + Ax[2][2]*Ax[2][2];
		Axx34 = Ax[2][1]*Ax[1][3];
		Axx41 = Ax[3][1]*Ax[1][0] + Ax[3][2]*Ax[2][0] + Ax[3][3]*Ax[3][0];
		Axx42 =      Ax[3][0]     + Ax[3][1]*Ax[1][1] + Ax[3][2]*Ax[2][1] + Ax[3][3]*Ax[3][1];
		Axx43 = Ax[3][1]*Ax[1][2] + Ax[3][2]*Ax[2][2] + Ax[3][3]*Ax[3][2];
		Axx44 = Ax[3][1]*Ax[1][3] + Ax[3][3]*Ax[3][3];

		double Axy11, Axy12, Axy13, Axy14, Axy21, Axy22, Axy23, Axy24, Axy31, Axy32, Axy33, Axy34, Axy41, Axy42, Axy43, Axy44;
		// calcula Axy = AxAy
		Axy11 = Ay[1][0];
		Axy12 = Ay[1][1];
		Axy13 = Ay[1][2];
		Axy14 = 0.0;
		Axy21 = Ax[1][1]*Ay[1][0] + Ax[1][2]*Ay[2][0] + Ax[1][3]*Ay[3][0];
		Axy22 = Ax[1][1]*Ay[1][1] + Ax[1][2]*Ay[2][1] + Ax[1][3]*Ay[3][1];
		Axy23 =      Ax[1][0]     + Ax[1][1]*Ay[1][2] + Ax[1][2]*Ay[2][2] + Ax[1][3]*Ay[3][2];
		Axy24 = Ax[1][2]*Ay[2][3] + Ax[1][3]*Ay[3][3];
		Axy31 = Ax[2][1]*Ay[1][0] + Ax[2][2]*Ay[2][0];
		Axy32 = Ax[2][1]*Ay[1][1] + Ax[2][2]*Ay[2][1];
		Axy33 =      Ax[2][0]     + Ax[2][1]*Ay[1][2] + Ax[2][2]*Ay[2][2];
		Axy34 = Ax[2][2]*Ay[2][3];
		Axy41 = Ax[3][1]*Ay[1][0] + Ax[3][2]*Ay[2][0] + Ax[3][3]*Ay[3][0];
		Axy42 = Ax[3][1]*Ay[1][1] + Ax[3][2]*Ay[2][1] + Ax[3][3]*Ay[3][1];
		Axy43 =      Ax[3][0]     + Ax[3][1]*Ay[1][2] + Ax[3][2]*Ay[2][2] + Ax[3][3]*Ay[3][2];
		Axy44 = Ax[3][2]*Ay[2][3] + Ax[3][3]*Ay[3][3];

		double Ayx11, Ayx12, Ayx13, Ayx14, Ayx21, Ayx22, Ayx23, Ayx24, Ayx31, Ayx32, Ayx33, Ayx34, Ayx41, Ayx42, Ayx43, Ayx44;
		//calcula Ayx = AyAx
		Ayx11 = Ax[2][0];
		Ayx12 = Ax[2][1];
		Ayx13 = Ax[2][2];
		Ayx14 = 0.0;
		Ayx21 = Ay[1][1]*Ax[1][0] + Ay[1][2]*Ax[2][0];
		Ayx22 =      Ay[1][0]     + Ay[1][1]*Ax[1][1] + Ay[1][2]*Ax[2][1];
		Ayx23 = Ay[1][1]*Ax[1][2] + Ay[1][2]*Ax[2][2];
		Ayx24 = Ay[1][1]*Ax[1][3];
		Ayx31 = Ay[2][1]*Ax[1][0] + Ay[2][2]*Ax[2][0] + Ay[2][3]*Ax[3][0];
		Ayx32 =      Ay[2][0]     + Ay[2][1]*Ax[1][1] + Ay[2][2]*Ax[2][1] + Ay[2][3]*Ax[3][1];
		Ayx33 = Ay[2][1]*Ax[1][2] + Ay[2][2]*Ax[2][2] + Ay[2][3]*Ax[3][2];
		Ayx34 = Ay[2][1]*Ax[1][3] + Ay[2][3]*Ax[3][3];
		Ayx41 = Ay[3][1]*Ax[1][0] + Ay[3][2]*Ax[2][0] + Ay[3][3]*Ax[3][0];
		Ayx42 =      Ay[3][0]     + Ay[3][1]*Ax[1][1] + Ay[3][2]*Ax[2][1] + Ay[3][3]*Ax[3][1];
		Ayx43 = Ay[3][1]*Ax[1][2] + Ay[3][2]*Ax[2][2] + Ay[3][3]*Ax[3][2];
		Ayx44 = Ay[3][1]*Ax[1][3] + Ay[3][3]*Ax[3][3];

		double Ayy11, Ayy12, Ayy13, Ayy14, Ayy21, Ayy22, Ayy23, Ayy24, Ayy31, Ayy32, Ayy33, Ayy34, Ayy41, Ayy42, Ayy43, Ayy44;
		//calcula Ayy = AyAy
		Ayy11 = Ay[2][0];
		Ayy12 = Ay[2][1];
		Ayy13 = Ay[2][2];
		Ayy14 = Ay[2][3];
		Ayy21 = Ay[1][1]*Ay[1][0] + Ay[1][2]*Ay[2][0];
		Ayy22 = Ay[1][1]*Ay[1][1] + Ay[1][2]*Ay[2][1];
		Ayy23 =      Ay[1][0]     + Ay[1][1]*Ay[1][2] + Ay[1][2]*Ay[2][2];
		Ayy24 = Ay[1][2]*Ay[2][3];
		Ayy31 = Ay[2][1]*Ay[1][0] + Ay[2][2]*Ay[2][0] + Ay[2][3]*Ay[3][0];
		Ayy32 = Ay[2][1]*Ay[1][1] + Ay[2][2]*Ay[2][1] + Ay[2][3]*Ay[3][1];
		Ayy33 =      Ay[2][0]     + Ay[2][1]*Ay[1][2] + Ay[2][2]*Ay[2][2] + Ay[2][3]*Ay[3][2];
		Ayy34 = Ay[2][2]*Ay[2][3] + Ay[2][3]*Ay[3][3];
		Ayy41 = Ay[3][1]*Ay[1][0] + Ay[3][2]*Ay[2][0] + Ay[3][3]*Ay[3][0];
		Ayy42 = Ay[3][1]*Ay[1][1] + Ay[3][2]*Ay[2][1] + Ay[3][3]*Ay[3][1];
		Ayy43 =      Ay[3][0]     + Ay[3][1]*Ay[1][2] + Ay[3][2]*Ay[2][2] + Ay[3][3]*Ay[3][2];
		Ayy44 = Ay[3][2]*Ay[2][3] + Ay[3][3]*Ay[3][3];

		//  OPERADOR DE CAPTURA DE DESCONTINUIDADES
		delta = FemFunctions->ShockCapture(tolerance, delta_old, gradUx, gradUy, Ax, Ay, A0, dUb, y23, y31, y12, x32, x13, x21,
											twoArea, e, Parameters->invY, Ub);

		// *** Calculo do parametro de estabilizacao (tau) do SUPG            
		betax = 2.0 * (gradUx[0]*Ub[0] + gradUx[1]*Ub[1] + gradUx[2]*Ub[2] + gradUx[3]*Ub[3]);
		betay = 2.0 * (gradUy[0]*Ub[0] + gradUy[1]*Ub[1] + gradUy[2]*Ub[2] + gradUy[3]*Ub[3]);

		betaxy = sqrt(betax*betax + betay*betay);    

		if (betaxy != 0.0)
		{
			betax = betax / betaxy;
			betay = betay / betaxy;
		}

		h = sqrt (twoArea);
		abs_vbeta = fabs(v1*betax + v2*betay);
		c = sqrt(gamma*(gamma-1.0)*(E - 0.5*(pow(v1,2) + pow(v2,2))));

		CFL = ((c + abs_vbeta) * delta_t) / h;
		psi = (2.0 * alpha * CFL) * (1.0/(1.0 + 2.0 * alpha * CFL));

		tau_d = delta / pow((c + abs_vbeta),2);
		tau_a = h / (2.0 * (c + abs_vbeta));
		tau_t = (2.0 / (3.0 * (1.0 + 2.0 * alpha * CFL))) * tau_a;

		tau = tau_t + psi * (tau_a - tau_d); 

		if (tau < 0.0) 
			tau = 0.0;

		//------------------------------------------------------------------------------
		//  MONTAGENS DAS MATRIZES
		//------------------------------------------------------------------------------

		// Coeficientes utilizados na matriz Mpg e Kg

		// Matriz Ca = y23Ax + x32Ay
		double Ca11, Ca12, Ca13, Ca14, Ca21, Ca22, Ca23, Ca24, Ca31, Ca32, Ca33, Ca34, Ca41, Ca42, Ca43, Ca44;
		Ca11 = 0.0;
		Ca12 = y23 * sixth;
		Ca13 = x32 * sixth;
		Ca14 = 0.0;
		Ca21 = (y23 * Ax[1][0] + x32 * Ay[1][0]) * sixth;
		Ca22 = (y23 * Ax[1][1] + x32 * Ay[1][1]) * sixth;
		Ca23 = (y23 * Ax[1][2] + x32 * Ay[1][2]) * sixth;
		Ca24 = y23 * Ax[1][3] * sixth;
		Ca31 = (y23 * Ax[2][0] + x32 * Ay[2][0]) * sixth;
		Ca32 = (y23 * Ax[2][1] + x32 * Ay[2][1]) * sixth;
		Ca33 = (y23 * Ax[2][2] + x32 * Ay[2][2]) * sixth;
		Ca34 = x32 * Ay[2][3] * sixth;
		Ca41 = (y23 * Ax[3][0] + x32 * Ay[3][0]) * sixth;
		Ca42 = (y23 * Ax[3][1] + x32 * Ay[3][1]) * sixth;
		Ca43 = (y23 * Ax[3][2] + x32 * Ay[3][2]) * sixth;
		Ca44 = (y23 * Ax[3][3] + x32 * Ay[3][3]) * sixth;

		// Matriz Cb = y31Ax + x13Ay
		double Cb11, Cb12, Cb13, Cb14, Cb21, Cb22, Cb23, Cb24, Cb31, Cb32, Cb33, Cb34, Cb41, Cb42, Cb43, Cb44;
		Cb11 = 0.0;
		Cb12 = y31 * sixth;
		Cb13 = x13 * sixth;
		Cb14 = 0.0;
		Cb21 = (y31 * Ax[1][0] + x13 * Ay[1][0]) * sixth;
		Cb22 = (y31 * Ax[1][1] + x13 * Ay[1][1]) * sixth;
		Cb23 = (y31 * Ax[1][2] + x13 * Ay[1][2]) * sixth;
		Cb24 = y31 * Ax[1][3] * sixth;
		Cb31 = (y31 * Ax[2][0] + x13 * Ay[2][0]) * sixth;
		Cb32 = (y31 * Ax[2][1] + x13 * Ay[2][1]) * sixth;
		Cb33 = (y31 * Ax[2][2] + x13 * Ay[2][2]) * sixth;
		Cb34 = x13 * Ay[2][3] * sixth;
		Cb41 = (y31 * Ax[3][0] + x13 * Ay[3][0]) * sixth;
		Cb42 = (y31 * Ax[3][1] + x13 * Ay[3][1]) * sixth;
		Cb43 = (y31 * Ax[3][2] + x13 * Ay[3][2]) * sixth;
		Cb44 = (y31 * Ax[3][3] + x13 * Ay[3][3]) * sixth;

		// Matriz Cc = y12Ax + x21Ay
		double Cc11, Cc12, Cc13, Cc14, Cc21, Cc22, Cc23, Cc24, Cc31, Cc32, Cc33, Cc34, Cc41, Cc42, Cc43, Cc44;
		Cc11 = 0.0;
		Cc12 = y12 * sixth;
		Cc13 = x21 * sixth;
		Cc14 = 0.0;
		Cc21 = (y12 * Ax[1][0] + x21 * Ay[1][0]) * sixth;
		Cc22 = (y12 * Ax[1][1] + x21 * Ay[1][1]) * sixth;
		Cc23 = (y12 * Ax[1][2] + x21 * Ay[1][2]) * sixth;
		Cc24 = y12 * Ax[1][3] * sixth;
		Cc31 = (y12 * Ax[2][0] + x21 * Ay[2][0]) * sixth;
		Cc32 = (y12 * Ax[2][1] + x21 * Ay[2][1]) * sixth;
		Cc33 = (y12 * Ax[2][2] + x21 * Ay[2][2]) * sixth;
		Cc34 = x21 * Ay[2][3] * sixth;
		Cc41 = (y12 * Ax[3][0] + x21 * Ay[3][0]) * sixth;
		Cc42 = (y12 * Ax[3][1] + x21 * Ay[3][1]) * sixth;
		Cc43 = (y12 * Ax[3][2] + x21 * Ay[3][2]) * sixth;
		Cc44 = (y12 * Ax[3][3] + x21 * Ay[3][3]) * sixth;

		// *** Matriz de Massa do Galerkin     
		double Mg1, Mg2;
		Mg1 = Area * sixth;
		Mg2 = Area / 12.0;

		// *** Matriz de Massa do SUPG 
		double Mpg11, Mpg12, Mpg13, Mpg14, Mpg21, Mpg22, Mpg23, Mpg24, Mpg31, Mpg32, Mpg33, Mpg34, Mpg41, Mpg42, Mpg43, Mpg44;
		double Mpg51, Mpg52, Mpg53, Mpg54, Mpg61, Mpg62, Mpg63, Mpg64, Mpg71, Mpg72, Mpg73, Mpg74, Mpg81, Mpg82, Mpg83, Mpg84;
		double Mpg91, Mpg92, Mpg93, Mpg94, Mpg101, Mpg102, Mpg103, Mpg104, Mpg111, Mpg112, Mpg113, Mpg114, Mpg121, Mpg122, Mpg123, Mpg124;

		// Matriz y23Ax + x32Ay
		Mpg11 = Ca11 * tau;
		Mpg12 = Ca12 * tau;
		Mpg13 = Ca13 * tau;
		Mpg14 = Ca14 * tau;
		Mpg21 = Ca21 * tau;
		Mpg22 = Ca22 * tau;
		Mpg23 = Ca23 * tau;
		Mpg24 = Ca24 * tau;
		Mpg31 = Ca31 * tau;
		Mpg32 = Ca32 * tau;
		Mpg33 = Ca33 * tau;
		Mpg34 = Ca34 * tau;
		Mpg41 = Ca41 * tau;
		Mpg42 = Ca42 * tau;
		Mpg43 = Ca43 * tau;
		Mpg44 = Ca44 * tau;

		// Matriz y31Ax + x13Ay
		Mpg51 = Cb11 * tau;
		Mpg52 = Cb12 * tau;
		Mpg53 = Cb13 * tau;
		Mpg54 = Cb14 * tau;
		Mpg61 = Cb21 * tau;
		Mpg62 = Cb22 * tau;
		Mpg63 = Cb23 * tau;
		Mpg64 = Cb24 * tau;
		Mpg71 = Cb31 * tau;
		Mpg72 = Cb32 * tau;
		Mpg73 = Cb33 * tau;
		Mpg74 = Cb34 * tau;
		Mpg81 = Cb41 * tau;
		Mpg82 = Cb42 * tau;
		Mpg83 = Cb43 * tau;
		Mpg84 = Cb44 * tau;

		// Matriz y12Ax + x21Ay
		Mpg91 = Cc11 * tau;
		Mpg92 = Cc12 * tau;
		Mpg93 = Cc13 * tau;
		Mpg94 = Cc14 * tau;
		Mpg101 = Cc21 * tau;
		Mpg102 = Cc22 * tau;
		Mpg103 = Cc23 * tau;
		Mpg104 = Cc24 * tau;
		Mpg111 = Cc31 * tau;
		Mpg112 = Cc32 * tau;
		Mpg113 = Cc33 * tau;
		Mpg114 = Cc34 * tau;
		Mpg121 = Cc41 * tau;
		Mpg122 = Cc42 * tau;
		Mpg123 = Cc43 * tau;
		Mpg124 = Cc44 * tau;

		 
		// *** Coeficientes da Matriz de Massa [Me]

		//no1 - no1
		Me[0][0] = Mg1 + Mpg11;
		Me[0][1] =       Mpg12;
		Me[0][2] =       Mpg13;
		Me[0][3] =       Mpg14;
		Me[1][0] =       Mpg21;
		Me[1][1] = Mg1 + Mpg22;
		Me[1][2] =       Mpg23;
		Me[1][3] =       Mpg24;
		Me[2][0] =       Mpg31;
		Me[2][1] =       Mpg32;
		Me[2][2] = Mg1 + Mpg33;
		Me[2][3] =       Mpg34;
		Me[3][0] =       Mpg41;
		Me[3][1] =       Mpg42;
		Me[3][2] =       Mpg43;
		Me[3][3] = Mg1 + Mpg44;

		//no1 - no2 
		Me[0][4] = Mg2 + Mpg11;
		Me[0][5] =       Mpg12;
		Me[0][6] =       Mpg13;
		Me[0][7] =       Mpg14;
		Me[1][4] =       Mpg21;
		Me[1][5] = Mg2 + Mpg22;
		Me[1][6] =       Mpg23;
		Me[1][7] =       Mpg24;
		Me[2][4] =       Mpg31;
		Me[2][5] =       Mpg32;
		Me[2][6] = Mg2 + Mpg33;
		Me[2][7] =       Mpg34;
		Me[3][4] =       Mpg41;
		Me[3][5] =       Mpg42;
		Me[3][6] =       Mpg43;
		Me[3][7] = Mg2 + Mpg44;

		//no1 - no3
		Me[0][8] = Mg2 + Mpg11;
		Me[0][9] =       Mpg12;
		Me[0][10] =       Mpg13;
		Me[0][11] =       Mpg14;
		Me[1][8] =       Mpg21;
		Me[1][9] = Mg2 + Mpg22;
		Me[1][10] =       Mpg23;
		Me[1][11] =       Mpg24;
		Me[2][8] =       Mpg31;
		Me[2][9] =       Mpg32;
		Me[2][10] = Mg2 + Mpg33;
		Me[2][11] =       Mpg34;
		Me[3][8] =       Mpg41;
		Me[3][9] =       Mpg42;
		Me[3][10] =       Mpg43;
		Me[3][11] = Mg2 + Mpg44;

		//no2 - no1
		Me[4][0] = Mg2 + Mpg51;
		Me[4][1] =       Mpg52;
		Me[4][2] =       Mpg53;
		Me[4][3] =       Mpg54;
		Me[5][0] =       Mpg61;
		Me[5][1] = Mg2 + Mpg62;
		Me[5][2] =       Mpg63;
		Me[5][3] =       Mpg64;
		Me[6][0] =       Mpg71;
		Me[6][1] =       Mpg72;
		Me[6][2] = Mg2 + Mpg73;
		Me[6][3] =       Mpg74;
		Me[7][0] =       Mpg81;
		Me[7][1] =       Mpg82;
		Me[7][2] =       Mpg83;
		Me[7][3] = Mg2 + Mpg84;

		//no2 - no2
		Me[4][4] = Mg1 + Mpg51;
		Me[4][5] =       Mpg52;
		Me[4][6] =       Mpg53;
		Me[4][7] =       Mpg54;
		Me[5][4] =       Mpg61;
		Me[5][5] = Mg1 + Mpg62;
		Me[5][6] =       Mpg63;
		Me[5][7] =       Mpg64;
		Me[6][4] =       Mpg71;
		Me[6][5] =       Mpg72;
		Me[6][6] = Mg1 + Mpg73;
		Me[6][7] =       Mpg74;
		Me[7][4] =       Mpg81;
		Me[7][5] =       Mpg82;
		Me[7][6] =       Mpg83;
		Me[7][7] = Mg1 + Mpg84;

		//no2 - no3
		Me[4][8] = Mg2 + Mpg51;
		Me[4][9] =       Mpg52;
		Me[4][10] =       Mpg53;
		Me[4][11] =       Mpg54;
		Me[5][8] =       Mpg61;
		Me[5][9] = Mg2 + Mpg62;
		Me[5][10] =       Mpg63;
		Me[5][11] =       Mpg64;
		Me[6][8] =       Mpg71;
		Me[6][9] =       Mpg72;
		Me[6][10] = Mg2 + Mpg73;
		Me[6][11] =       Mpg74;
		Me[7][8] =       Mpg81;
		Me[7][9] =       Mpg82;
		Me[7][10] =       Mpg83;
		Me[7][11] = Mg2 + Mpg84;

		//no3 - no1
		Me[8][0] = Mg2 + Mpg91;
		Me[8][1] =       Mpg92;
		Me[8][2] =       Mpg93;
		Me[8][3] =       Mpg94;
		Me[9][0] =       Mpg101;
		Me[9][1] = Mg2 + Mpg102;
		Me[9][2] =       Mpg103;
		Me[9][3] =       Mpg104;
		Me[10][0] =       Mpg111;
		Me[10][1] =       Mpg112;
		Me[10][2] = Mg2 + Mpg113;
		Me[10][3] =       Mpg114;
		Me[11][0] =       Mpg121;
		Me[11][1] =       Mpg122;
		Me[11][2] =       Mpg123;
		Me[11][3] = Mg2 + Mpg124;

		//no3 - no2
		Me[8][4] = Mg2 + Mpg91;
		Me[8][5] =       Mpg92;
		Me[8][6] =       Mpg93;
		Me[8][7] =       Mpg94;
		Me[9][4] =       Mpg101;
		Me[9][5] = Mg2 + Mpg102;
		Me[9][6] =       Mpg103;
		Me[9][7] =       Mpg104;
		Me[10][4] =       Mpg111;
		Me[10][5] =       Mpg112;
		Me[10][6] = Mg2 + Mpg113;
		Me[10][7] =       Mpg114;
		Me[11][4] =       Mpg121;
		Me[11][5] =       Mpg122;
		Me[11][6] =       Mpg123;
		Me[11][7] = Mg2 + Mpg124;

		//no3 - no3
		Me[8][8] = Mg1 + Mpg91;
		Me[8][9] =       Mpg92;
		Me[8][10] =       Mpg93;
		Me[8][11] =       Mpg94;
		Me[9][8] =       Mpg101;
		Me[9][9] = Mg1 + Mpg102;
		Me[9][10] =       Mpg103;
		Me[9][11] =       Mpg104;
		Me[10][8] =       Mpg111;
		Me[10][9] =       Mpg112;
		Me[10][10] = Mg1 + Mpg113;
		Me[10][11] =       Mpg114;
		Me[11][8] =       Mpg121;
		Me[11][9] =       Mpg122;
		Me[11][10] =       Mpg123;
		Me[11][11] = Mg1 + Mpg124;

		// *** Matriz de Conveccao de Galerkin
		double Kg11, Kg12, Kg13, Kg14, Kg21, Kg22, Kg23, Kg24, Kg31, Kg32, Kg33, Kg34, Kg41, Kg42, Kg43, Kg44;
		double Kg15, Kg16, Kg17, Kg18, Kg25, Kg26, Kg27, Kg28, Kg35, Kg36, Kg37, Kg38, Kg45, Kg46, Kg47, Kg48;
		double Kg19, Kg110, Kg111, Kg112, Kg29, Kg210, Kg211, Kg212, Kg39, Kg310, Kg311, Kg312, Kg49, Kg410, Kg411, Kg412;

		Kg11 = Ca11;
		Kg12 = Ca12;
		Kg13 = Ca13;
		Kg14 = Ca14;
		Kg21 = Ca21;
		Kg22 = Ca22;
		Kg23 = Ca23;
		Kg24 = Ca24;
		Kg31 = Ca31;
		Kg32 = Ca32;
		Kg33 = Ca33;
		Kg34 = Ca34;
		Kg41 = Ca41;
		Kg42 = Ca42;
		Kg43 = Ca43;
		Kg44 = Ca44;

		Kg15 = Cb11;
		Kg16 = Cb12;
		Kg17 = Cb13;
		Kg18 = Cb14;
		Kg25 = Cb21;
		Kg26 = Cb22;
		Kg27 = Cb23;
		Kg28 = Cb24;
		Kg35 = Cb31;
		Kg36 = Cb32;
		Kg37 = Cb33;
		Kg38 = Cb34;
		Kg45 = Cb41;
		Kg46= Cb42;
		Kg47 = Cb43;
		Kg48 = Cb44;

		Kg19 = Cc11;
		Kg110 = Cc12;
		Kg111 = Cc13;
		Kg112 = Cc14;
		Kg29 = Cc21;
		Kg210 = Cc22;
		Kg211 = Cc23;
		Kg212 = Cc24;
		Kg39 = Cc31;
		Kg310 = Cc32;
		Kg311 = Cc33;
		Kg312 = Cc34;
		Kg49 = Cc41;
		Kg410 = Cc42;
		Kg411 = Cc43;
		Kg412 = Cc44;
		  
		// *** Matriz de Correcao SUPG 

		//  coeficientes da matriz resultante da multiplicação B^t[AxAx  AxAy // AyAx  AyAy]
		double a23x11, a23x12, a23x13, a23x14, a23x21, a23x22, a23x23, a23x24, a23x31, a23x32, a23x33, a23x34, a23x41, a23x42, a23x43, a23x44;
		double a23y11, a23y12, a23y13, a23y14, a23y21, a23y22, a23y23, a23y24, a23y31, a23y32, a23y33, a23y34, a23y41, a23y42, a23y43, a23y44;
		double a31x11, a31x12, a31x13, a31x14, a31x21, a31x22, a31x23, a31x24, a31x31, a31x32, a31x33, a31x34, a31x41, a31x42, a31x43, a31x44;
		double a31y11, a31y12, a31y13, a31y14, a31y21, a31y22, a31y23, a31y24, a31y31, a31y32, a31y33, a31y34, a31y41, a31y42, a31y43, a31y44;
		double a12x11, a12x12, a12x13, a12x14, a12x21, a12x22, a12x23, a12x24, a12x31, a12x32, a12x33, a12x34, a12x41, a12x42, a12x43, a12x44;
		double a12y11, a12y12, a12y13, a12y14, a12y21, a12y22, a12y23, a12y24, a12y31, a12y32, a12y33, a12y34, a12y41, a12y42, a12y43, a12y44;

		// y23*Axx + x32*Ayx = [a23xij]
		a23x11 = y23*Axx11 + x32*Ayx11;
		a23x12 = y23*Axx12 + x32*Ayx12;
		a23x13 = y23*Axx13 + x32*Ayx13;
		a23x14 = y23*Axx14 + x32*Ayx14;
		a23x21 = y23*Axx21 + x32*Ayx21;
		a23x22 = y23*Axx22 + x32*Ayx22;
		a23x23 = y23*Axx23 + x32*Ayx23;
		a23x24 = y23*Axx24 + x32*Ayx24;
		a23x31 = y23*Axx31 + x32*Ayx31;
		a23x32 = y23*Axx32 + x32*Ayx32;
		a23x33 = y23*Axx33 + x32*Ayx33;
		a23x34 = y23*Axx34 + x32*Ayx34;
		a23x41 = y23*Axx41 + x32*Ayx41;
		a23x42 = y23*Axx42 + x32*Ayx42;
		a23x43 = y23*Axx43 + x32*Ayx43;
		a23x44 = y23*Axx44 + x32*Ayx44;
                                  
		// y23*Axy + x32*Ayy = [a23yij]
		a23y11 = y23*Axy11 + x32*Ayy11;
		a23y12 = y23*Axy12 + x32*Ayy12;
		a23y13 = y23*Axy13 + x32*Ayy13;
		a23y14 = y23*Axy14 + x32*Ayy14;
		a23y21 = y23*Axy21 + x32*Ayy21;
		a23y22 = y23*Axy22 + x32*Ayy22;
		a23y23 = y23*Axy23 + x32*Ayy23;
		a23y24 = y23*Axy24 + x32*Ayy24;
		a23y31 = y23*Axy31 + x32*Ayy31;
		a23y32 = y23*Axy32 + x32*Ayy32;
		a23y33 = y23*Axy33 + x32*Ayy33;
		a23y34 = y23*Axy34 + x32*Ayy34;
		a23y41 = y23*Axy41 + x32*Ayy41;
		a23y42 = y23*Axy42 + x32*Ayy42;
		a23y43 = y23*Axy43 + x32*Ayy43;
		a23y44 = y23*Axy44 + x32*Ayy44;

		// y31*Axx + x13*Ayx = [a31xij]
		a31x11 = y31*Axx11 + x13*Ayx11;
		a31x12 = y31*Axx12 + x13*Ayx12;
		a31x13 = y31*Axx13 + x13*Ayx13;
		a31x14 = y31*Axx14 + x13*Ayx14;
		a31x21 = y31*Axx21 + x13*Ayx21;
		a31x22 = y31*Axx22 + x13*Ayx22;
		a31x23 = y31*Axx23 + x13*Ayx23;
		a31x24 = y31*Axx24 + x13*Ayx24;
		a31x31 = y31*Axx31 + x13*Ayx31;
		a31x32 = y31*Axx32 + x13*Ayx32;
		a31x33 = y31*Axx33 + x13*Ayx33;
		a31x34 = y31*Axx34 + x13*Ayx34;
		a31x41 = y31*Axx41 + x13*Ayx41;
		a31x42 = y31*Axx42 + x13*Ayx42;
		a31x43 = y31*Axx43 + x13*Ayx43;
		a31x44 = y31*Axx44 + x13*Ayx44;

		// y31*Axy + x13*Ayy = [a31yij]
		a31y11 = y31*Axy11 + x13*Ayy11;
		a31y12 = y31*Axy12 + x13*Ayy12;
		a31y13 = y31*Axy13 + x13*Ayy13;    
		a31y14 = y31*Axy14 + x13*Ayy14;
		a31y21 = y31*Axy21 + x13*Ayy21;
		a31y22 = y31*Axy22 + x13*Ayy22;
		a31y23 = y31*Axy23 + x13*Ayy23;
		a31y24 = y31*Axy24 + x13*Ayy24;
		a31y31 = y31*Axy31 + x13*Ayy31;
		a31y32 = y31*Axy32 + x13*Ayy32;
		a31y33 = y31*Axy33 + x13*Ayy33;
		a31y34 = y31*Axy34 + x13*Ayy34;
		a31y41 = y31*Axy41 + x13*Ayy41;
		a31y42 = y31*Axy42 + x13*Ayy42;
		a31y43 = y31*Axy43 + x13*Ayy43;
		a31y44 = y31*Axy44 + x13*Ayy44;

		// y12*Axx + x21*Ayx = [a12xij]
		a12x11 = y12*Axx11 + x21*Ayx11;
		a12x12 = y12*Axx12 + x21*Ayx12;
		a12x13 = y12*Axx13 + x21*Ayx13;    
		a12x14 = y12*Axx14 + x21*Ayx14;
		a12x21 = y12*Axx21 + x21*Ayx21;
		a12x22 = y12*Axx22 + x21*Ayx22;
		a12x23 = y12*Axx23 + x21*Ayx23;
		a12x24 = y12*Axx24 + x21*Ayx24;
		a12x31 = y12*Axx31 + x21*Ayx31;
		a12x32 = y12*Axx32 + x21*Ayx32;
		a12x33 = y12*Axx33 + x21*Ayx33;
		a12x34 = y12*Axx34 + x21*Ayx34;
		a12x41 = y12*Axx41 + x21*Ayx41;
		a12x42 = y12*Axx42 + x21*Ayx42;
		a12x43 = y12*Axx43 + x21*Ayx43;
		a12x44 = y12*Axx44 + x21*Ayx44;

		// y12*Axy + x21*Ayy = [a12yij]
		a12y11 = y12*Axy11 + x21*Ayy11;
		a12y12 = y12*Axy12 + x21*Ayy12;
		a12y13 = y12*Axy13 + x21*Ayy13;    
		a12y14 = y12*Axy14 + x21*Ayy14;
		a12y21 = y12*Axy21 + x21*Ayy21;
		a12y22 = y12*Axy22 + x21*Ayy22;
		a12y23 = y12*Axy23 + x21*Ayy23;
		a12y24 = y12*Axy24 + x21*Ayy24;
		a12y31 = y12*Axy31 + x21*Ayy31;
		a12y32 = y12*Axy32 + x21*Ayy32;
		a12y33 = y12*Axy33 + x21*Ayy33;
		a12y34 = y12*Axy34 + x21*Ayy34;
		a12y41 = y12*Axy41 + x21*Ayy41;
		a12y42 = y12*Axy42 + x21*Ayy42;
		a12y43 = y12*Axy43 + x21*Ayy43;
		a12y44 = y12*Axy44 + x21*Ayy44;

		// Coeficiente para a montagem da Kpg - não precisa calcular os termos da diagonal
		double Kpg15, Kpg16, Kpg17, Kpg18, Kpg25, Kpg26, Kpg27, Kpg28, Kpg35, Kpg36, Kpg37, Kpg38, Kpg45, Kpg46, Kpg47, Kpg48;
		double Kpg19, Kpg110, Kpg1_11, Kpg1_12, Kpg29, Kpg210, Kpg211, Kpg212, Kpg39, Kpg310, Kpg311, Kpg312, Kpg49, Kpg410, Kpg411, Kpg412;
		double Kpg51, Kpg52, Kpg53, Kpg54, Kpg61, Kpg62, Kpg63, Kpg64, Kpg71, Kpg72, Kpg73, Kpg74, Kpg81, Kpg82, Kpg83, Kpg84;
		double Kpg59, Kpg510, Kpg511, Kpg512, Kpg69, Kpg610, Kpg611, Kpg612, Kpg79, Kpg710, Kpg711, Kpg712, Kpg89, Kpg810, Kpg811, Kpg812;
		double Kpg91, Kpg92, Kpg93, Kpg94, Kpg101, Kpg102, Kpg103, Kpg104, Kpg11_1, Kpg11_2, Kpg113, Kpg114, Kpg121, Kpg122, Kpg123, Kpg124;
		double Kpg95, Kpg96, Kpg97, Kpg98, Kpg105, Kpg106, Kpg107, Kpg108, Kpg115, Kpg116, Kpg117, Kpg118, Kpg125, Kpg126, Kpg127, Kpg128;

		double tau4area = tau / (4 * Area);

		// y31*a23x + x13*a23y
		Kpg15 = tau4area * (y31*a23x11 + x13*a23y11);
		Kpg16 = tau4area * (y31*a23x12 + x13*a23y12);
		Kpg17 = tau4area * (y31*a23x13 + x13*a23y13);
		Kpg18 = tau4area * (y31*a23x14 + x13*a23y14);
		Kpg25 = tau4area * (y31*a23x21 + x13*a23y21);
		Kpg26 = tau4area * (y31*a23x22 + x13*a23y22);
		Kpg27 = tau4area * (y31*a23x23 + x13*a23y23);
		Kpg28 = tau4area * (y31*a23x24 + x13*a23y24);
		Kpg35 = tau4area * (y31*a23x31 + x13*a23y31);
		Kpg36 = tau4area * (y31*a23x32 + x13*a23y32);
		Kpg37 = tau4area * (y31*a23x33 + x13*a23y33);
		Kpg38 = tau4area * (y31*a23x34 + x13*a23y34);
		Kpg45 = tau4area * (y31*a23x41 + x13*a23y41);
		Kpg46 = tau4area * (y31*a23x42 + x13*a23y42);
		Kpg47 = tau4area * (y31*a23x43 + x13*a23y43);
		Kpg48 = tau4area * (y31*a23x44 + x13*a23y44);

		// y12*a23x + x21*a23y
		Kpg19   = tau4area * (y12*a23x11 + x21*a23y11);
		Kpg110  = tau4area * (y12*a23x12 + x21*a23y12);
		Kpg1_11 = tau4area * (y12*a23x13 + x21*a23y13);
		Kpg1_12 = tau4area * (y12*a23x14 + x21*a23y14);
		Kpg29   = tau4area * (y12*a23x21 + x21*a23y21);
		Kpg210  = tau4area * (y12*a23x22 + x21*a23y22);
		Kpg211  = tau4area * (y12*a23x23 + x21*a23y23);
		Kpg212  = tau4area * (y12*a23x24 + x21*a23y24);
		Kpg39   = tau4area * (y12*a23x31 + x21*a23y31);
		Kpg310  = tau4area * (y12*a23x32 + x21*a23y32);
		Kpg311  = tau4area * (y12*a23x33 + x21*a23y33);
		Kpg312  = tau4area * (y12*a23x34 + x21*a23y34);
		Kpg49   = tau4area * (y12*a23x41 + x21*a23y41);
		Kpg410  = tau4area * (y12*a23x42 + x21*a23y42);
		Kpg411  = tau4area * (y12*a23x43 + x21*a23y43);
		Kpg412  = tau4area * (y12*a23x44 + x21*a23y44);

		// y23*a31x + x32*a31y
		Kpg51 = tau4area * (y23*a31x11 + x32*a31y11);
		Kpg52 = tau4area * (y23*a31x12 + x32*a31y12);
		Kpg53 = tau4area * (y23*a31x13 + x32*a31y13);
		Kpg54 = tau4area * (y23*a31x14 + x32*a31y14);
		Kpg61 = tau4area * (y23*a31x21 + x32*a31y21);
		Kpg62 = tau4area * (y23*a31x22 + x32*a31y22);
		Kpg63 = tau4area * (y23*a31x23 + x32*a31y23);
		Kpg64 = tau4area * (y23*a31x24 + x32*a31y24);
		Kpg71 = tau4area * (y23*a31x31 + x32*a31y31);
		Kpg72 = tau4area * (y23*a31x32 + x32*a31y32);
		Kpg73 = tau4area * (y23*a31x33 + x32*a31y33);
		Kpg74 = tau4area * (y23*a31x34 + x32*a31y34);
		Kpg81 = tau4area * (y23*a31x41 + x32*a31y41);
		Kpg82 = tau4area * (y23*a31x42 + x32*a31y42);
		Kpg83 = tau4area * (y23*a31x43 + x32*a31y43);
		Kpg84 = tau4area * (y23*a31x44 + x32*a31y44);

		// y12*a31x + x21*a31y
		Kpg59  = tau4area * (y12*a31x11 + x21*a31y11);
		Kpg510 = tau4area * (y12*a31x12 + x21*a31y12);
		Kpg511 = tau4area * (y12*a31x13 + x21*a31y13);
		Kpg512 = tau4area * (y12*a31x14 + x21*a31y14);
		Kpg69  = tau4area * (y12*a31x21 + x21*a31y21);
		Kpg610 = tau4area * (y12*a31x22 + x21*a31y22);
		Kpg611 = tau4area * (y12*a31x23 + x21*a31y23);
		Kpg612 = tau4area * (y12*a31x24 + x21*a31y24);
		Kpg79  = tau4area * (y12*a31x31 + x21*a31y31);
		Kpg710 = tau4area * (y12*a31x32 + x21*a31y32);
		Kpg711 = tau4area * (y12*a31x33 + x21*a31y33);
		Kpg712 = tau4area * (y12*a31x34 + x21*a31y34);
		Kpg89  = tau4area * (y12*a31x41 + x21*a31y41);
		Kpg810 = tau4area * (y12*a31x42 + x21*a31y42);
		Kpg811 = tau4area * (y12*a31x43 + x21*a31y43);
		Kpg812 = tau4area * (y12*a31x44 + x21*a31y44);

		// y23*a12x + x32*a12y
		Kpg91   = tau4area * (y23*a12x11 + x32*a12y11);
		Kpg92   = tau4area * (y23*a12x12 + x32*a12y12);
		Kpg93   = tau4area * (y23*a12x13 + x32*a12y13);
		Kpg94   = tau4area * (y23*a12x14 + x32*a12y14);
		Kpg101  = tau4area * (y23*a12x21 + x32*a12y21);
		Kpg102  = tau4area * (y23*a12x22 + x32*a12y22);
		Kpg103  = tau4area * (y23*a12x23 + x32*a12y23);
		Kpg104  = tau4area * (y23*a12x24 + x32*a12y24);
		Kpg11_1 = tau4area * (y23*a12x31 + x32*a12y31);
		Kpg11_2 = tau4area * (y23*a12x32 + x32*a12y32);
		Kpg113  = tau4area * (y23*a12x33 + x32*a12y33);
		Kpg114  = tau4area * (y23*a12x34 + x32*a12y34);
		Kpg121  = tau4area * (y23*a12x41 + x32*a12y41);
		Kpg122  = tau4area * (y23*a12x42 + x32*a12y42);
		Kpg123  = tau4area * (y23*a12x43 + x32*a12y43);
		Kpg124  = tau4area * (y23*a12x44 + x32*a12y44);

		// y31*a12x + x13*a12y
		Kpg95  = tau4area * (y31*a12x11 + x13*a12y11);
		Kpg96  = tau4area * (y31*a12x12 + x13*a12y12);
		Kpg97  = tau4area * (y31*a12x13 + x13*a12y13);
		Kpg98  = tau4area * (y31*a12x14 + x13*a12y14);
		Kpg105 = tau4area * (y31*a12x21 + x13*a12y21);
		Kpg106 = tau4area * (y31*a12x22 + x13*a12y22);
		Kpg107 = tau4area * (y31*a12x23 + x13*a12y23);
		Kpg108 = tau4area * (y31*a12x24 + x13*a12y24);
		Kpg115 = tau4area * (y31*a12x31 + x13*a12y31);
		Kpg116 = tau4area * (y31*a12x32 + x13*a12y32);
		Kpg117 = tau4area * (y31*a12x33 + x13*a12y33);
		Kpg118 = tau4area * (y31*a12x34 + x13*a12y34);
		Kpg125 = tau4area * (y31*a12x41 + x13*a12y41);
		Kpg126 = tau4area * (y31*a12x42 + x13*a12y42);
		Kpg127 = tau4area * (y31*a12x43 + x13*a12y43);
		Kpg128 = tau4area * (y31*a12x44 + x13*a12y44);
		 

		// *** Matriz de Correcao CAU
		double delta4area = delta / (4.0 * Area);

		double Kcau12 = delta4area * (y23*y31 + x32*x13);
		double Kcau13 = delta4area * (y23*y12 + x32*x21);
		double Kcau23 = delta4area * (y31*y12 + x13*x21);


		// *** Coeficientes da Matriz de Rigidez [Ke]

		//no1 - no1
		Ke[0][0] = Kg11 - (Kpg15 + Kpg19) - (Kcau12 + Kcau13);
		Ke[0][1] = Kg12 - (Kpg16 + Kpg110);
		Ke[0][2] = Kg13 - (Kpg17 + Kpg1_11);
		Ke[0][3] = Kg14 - (Kpg18 + Kpg1_12);
		Ke[1][0] = Kg21 - (Kpg25 + Kpg29);
		Ke[1][1] = Kg22 - (Kpg26 + Kpg210) - (Kcau12 + Kcau13);
		Ke[1][2] = Kg23 - (Kpg27 + Kpg211);
		Ke[1][3] = Kg24 - (Kpg28 + Kpg212);
		Ke[2][0] = Kg31 - (Kpg35 + Kpg39);
		Ke[2][1] = Kg32 - (Kpg36 + Kpg310);
		Ke[2][2] = Kg33 - (Kpg37 + Kpg311) - (Kcau12 + Kcau13);
		Ke[2][3] = Kg34 - (Kpg38 + Kpg312);
		Ke[3][0] = Kg41 - (Kpg45 + Kpg49);
		Ke[3][1] = Kg42 - (Kpg46 + Kpg410);
		Ke[3][2] = Kg43 - (Kpg47 + Kpg411);
		Ke[3][3] = Kg44 - (Kpg48 + Kpg412) - (Kcau12 + Kcau13);

		//no1 - no2
		Ke[0][4] = Kg15 + Kpg15 + Kcau12;
		Ke[0][5] = Kg16 + Kpg16;
		Ke[0][6] = Kg17 + Kpg17;
		Ke[0][7] = Kg18 + Kpg18;
		Ke[1][4] = Kg25 + Kpg25;
		Ke[1][5] = Kg26 + Kpg26 + Kcau12;
		Ke[1][6] = Kg27 + Kpg27;
		Ke[1][7] = Kg28 + Kpg28;
		Ke[2][4] = Kg35 + Kpg35;
		Ke[2][5] = Kg36 + Kpg36;
		Ke[2][6] = Kg37 + Kpg37 + Kcau12;
		Ke[2][7] = Kg38 + Kpg38;
		Ke[3][4] = Kg45 + Kpg45;
		Ke[3][5] = Kg46 + Kpg46;
		Ke[3][6] = Kg47 + Kpg47;
		Ke[3][7] = Kg48 + Kpg48 + Kcau12;

		//no1 - no3
		Ke[0][8]  = Kg19  + Kpg19 + Kcau13;
		Ke[0][9]  = Kg110 + Kpg110;
		Ke[0][10] = Kg111 + Kpg1_11;
		Ke[0][11] = Kg112 + Kpg1_12;
		Ke[1][8]  = Kg29  + Kpg29;
		Ke[1][9]  = Kg210 + Kpg210 + Kcau13;
		Ke[1][10] = Kg211 + Kpg211;
		Ke[1][11] = Kg212 + Kpg212;
		Ke[2][8]  = Kg39  + Kpg39;
		Ke[2][9]  = Kg310 + Kpg310;
		Ke[2][10] = Kg311 + Kpg311 + Kcau13;
		Ke[2][11] = Kg312 + Kpg312;
		Ke[3][8]  = Kg49  + Kpg49;
		Ke[3][9]  = Kg410 + Kpg410;
		Ke[3][10] = Kg411 + Kpg411;
		Ke[3][11] = Kg412 + Kpg412 + Kcau13;

		//no2 - no1
		Ke[4][0] = Kg11 + Kpg51 + Kcau12;
		Ke[4][1] = Kg12 + Kpg52;
		Ke[4][2] = Kg13 + Kpg53;
		Ke[4][3] = Kg14 + Kpg54;
		Ke[5][0] = Kg21 + Kpg61;
		Ke[5][1] = Kg22 + Kpg62 + Kcau12;
		Ke[5][2] = Kg23 + Kpg63;
		Ke[5][3] = Kg24 + Kpg64;
		Ke[6][0] = Kg31 + Kpg71;
		Ke[6][1] = Kg32 + Kpg72;
		Ke[6][2] = Kg33 + Kpg73 + Kcau12;
		Ke[6][3] = Kg34 + Kpg74;
		Ke[7][0] = Kg41 + Kpg81;
		Ke[7][1] = Kg42 + Kpg82;
		Ke[7][2] = Kg43 + Kpg83;
		Ke[7][3] = Kg44 + Kpg84 + Kcau12;

		//no2 - no2
		Ke[4][4] = Kg15 - (Kpg51 + Kpg59) - (Kcau12 + Kcau23);
		Ke[4][5] = Kg16 - (Kpg52 + Kpg510);
		Ke[4][6] = Kg17 - (Kpg53 + Kpg511);
		Ke[4][7] = Kg18 - (Kpg54 + Kpg512);
		Ke[5][4] = Kg25 - (Kpg61 + Kpg69);
		Ke[5][5] = Kg26 - (Kpg62 + Kpg610) - (Kcau12 + Kcau23);
		Ke[5][6] = Kg27 - (Kpg63 + Kpg611);
		Ke[5][7] = Kg28 - (Kpg64 + Kpg612);
		Ke[6][4] = Kg35 - (Kpg71 + Kpg79);
		Ke[6][5] = Kg36 - (Kpg72 + Kpg710);
		Ke[6][6] = Kg37 - (Kpg73 + Kpg711) - (Kcau12 + Kcau23);
		Ke[6][7] = Kg38 - (Kpg74 + Kpg712);
		Ke[7][4] = Kg45 - (Kpg81 + Kpg89);
		Ke[7][5] = Kg46 - (Kpg82 + Kpg810);
		Ke[7][6] = Kg47 - (Kpg83 + Kpg811);
		Ke[7][7] = Kg48 - (Kpg84 + Kpg812) - (Kcau12 + Kcau23);

		//no2 - no3
		Ke[4][8]  = Kg19  + Kpg59 + Kcau23;
		Ke[4][9]  = Kg110 + Kpg510;
		Ke[4][10] = Kg111 + Kpg511;
		Ke[4][11] = Kg112 + Kpg512;
		Ke[5][8]  = Kg29  + Kpg69;
		Ke[5][9]  = Kg210 + Kpg610 + Kcau23;
		Ke[5][10] = Kg211 + Kpg611;
		Ke[5][11] = Kg212 + Kpg612;
		Ke[6][8]  = Kg39  + Kpg79;
		Ke[6][9]  = Kg310 + Kpg710;
		Ke[6][10] = Kg311 + Kpg711 + Kcau23;
		Ke[6][11] = Kg312 + Kpg712;
		Ke[7][8]  = Kg49  + Kpg89;
		Ke[7][9]  = Kg410 + Kpg810;
		Ke[7][10] = Kg411 + Kpg811;
		Ke[7][11] = Kg412 + Kpg812 + Kcau23;

		//no3 - no1
		Ke[8][0]  = Kg11 + Kpg91 + Kcau13;
		Ke[8][1]  = Kg12 + Kpg92;
		Ke[8][2]  = Kg13 + Kpg93;
		Ke[8][3]  = Kg14 + Kpg94;
		Ke[9][0]  = Kg21 + Kpg101;
		Ke[9][1]  = Kg22 + Kpg102 + Kcau13;
		Ke[9][2]  = Kg23 + Kpg103;
		Ke[9][3]  = Kg24 + Kpg104;
		Ke[10][0] = Kg31 + Kpg11_1;
		Ke[10][1] = Kg32 + Kpg11_2;
		Ke[10][2] = Kg33 + Kpg113 + Kcau13;
		Ke[10][3] = Kg34 + Kpg114;
		Ke[11][0] = Kg41 + Kpg121;
		Ke[11][1] = Kg42 + Kpg122;
		Ke[11][2] = Kg43 + Kpg123;
		Ke[11][3] = Kg44 + Kpg124 + Kcau13;

		//no3 - no2
		Ke[8][4]  = Kg15 + Kpg95 + Kcau23;
		Ke[8][5]  = Kg16 + Kpg96;
		Ke[8][6]  = Kg17 + Kpg97;
		Ke[8][7]  = Kg18 + Kpg98;
		Ke[9][4]  = Kg25 + Kpg105;
		Ke[9][5]  = Kg26 + Kpg106 + Kcau23;
		Ke[9][6]  = Kg27 + Kpg107;
		Ke[9][7]  = Kg28 + Kpg108;
		Ke[10][4] = Kg35 + Kpg115;
		Ke[10][5] = Kg36 + Kpg116;
		Ke[10][6] = Kg37 + Kpg117 + Kcau23;
		Ke[10][7] = Kg38 + Kpg118;
		Ke[11][4] = Kg45 + Kpg125;
		Ke[11][5] = Kg46 + Kpg126;
		Ke[11][6] = Kg47 + Kpg127;
		Ke[11][7] = Kg48 + Kpg128 + Kcau23;

		//no3 - no3
		Ke[8][8]   = Kg19  - (Kpg91   + Kpg95) - (Kcau13 + Kcau23);
		Ke[8][9]   = Kg110 - (Kpg92   + Kpg96);
		Ke[8][10]  = Kg111 - (Kpg93   + Kpg97);
		Ke[8][11]  = Kg112 - (Kpg94   + Kpg98);
		Ke[9][8]   = Kg29  - (Kpg101  + Kpg105);
		Ke[9][9]   = Kg210 - (Kpg102  + Kpg106) - (Kcau13 + Kcau23);
		Ke[9][10]  = Kg211 - (Kpg103  + Kpg107);
		Ke[9][11]  = Kg212 - (Kpg104  + Kpg108);
		Ke[10][8]  = Kg39  - (Kpg11_1 + Kpg115);
		Ke[10][9]  = Kg310 - (Kpg11_2 + Kpg116);
		Ke[10][10] = Kg311 - (Kpg113  + Kpg117) - (Kcau13 + Kcau23);
		Ke[10][11] = Kg312 - (Kpg114  + Kpg118);
		Ke[11][8]  = Kg49  - (Kpg121  + Kpg125);
		Ke[11][9]  = Kg410 - (Kpg122  + Kpg126);
		Ke[11][10] = Kg411 - (Kpg123  + Kpg127);
		Ke[11][11] = Kg412 - (Kpg124  + Kpg128) - (Kcau13 + Kcau23); 

		//*****************************************

		//Fill local Ae and Re and Global R
		double Re[12], MedUe[12], KeUe[12], Ae[12][12];
		
		for (i=0; i<12; i++){
			MedUe[i] = 0;
			KeUe[i] = 0;
			for (j=0; j<12; j++){
				MedUe[i] += Me[i][j]*dUe[j];
				KeUe[i] += Ke[i][j]*Ue[j];
				Ae[i][j] = Me[i][j] + alpha*delta_t*Ke[i][j];
			}
			Re[i] = - MedUe[i] - KeUe[i];
			R[lm[e][i]] += Re[i];
		}
		
		R[NEQ] = 0;
		
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, e, Ae);

	}//for elemento

	// Liberacao dos espacos alocados
	free(U);
	free(dU);

	return 0;
}// end build

