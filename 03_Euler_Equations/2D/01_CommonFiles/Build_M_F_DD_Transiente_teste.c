#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_M_F_DD_Transiente(ParametersType *Parameters, double *delta_old, double *delta_old_NMV, double delta_t, double alpha, int **lm, MatrixDataType *MatrixData, double *F, double *u, double *du,
			 double *uB, double *duB, NodeType *Node, ElementType *Element, double (*ShockCapture)(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4], double *, double, double, double, double, double, double, double, int, double *, double *), double (*rhopresc)(double, double), double (*v1presc)(double, double), double (*v2presc)(double, double), double (*epresc)(double, double), double (*Gamma)(double, double), void (*assembly)(int, int, double (*)[12],  MatrixDataType *), AuxBuildStructuresType *AuxBuild, double t)
{
	int e, i, j, neq, nel, lm_aux[12];
	int J1, J2, J3;
	double x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21;
	double *vet_bound;
	double Area, twoArea, third = 1.0/3.0, sixth = 1.0/6.0, ninefortieth = 9.0 / 40.0;
	double *Ub, *Ue;
	double delta, gamma;
	double *gradUx, *gradUy;
	double *dUe, *dUb;
	double **M2_out, **R1_out, **R2_out, *invN2_out;
	double Me[12][12], Fe[12];
	double tolerance;// Mmin = 3.0, Mmax = 0.0;
	//char FileName[300];
	//FILE *outFile;
	
	tolerance = AuxBuild->tolerance;
	neq = Parameters->neq;
	nel = Parameters->nel;
	M2_out = AuxBuild->M2;
	R1_out = AuxBuild->R1;
	R2_out = AuxBuild->R2;
	invN2_out = AuxBuild->invN2;	
	
	vet_bound = (double*) mycalloc("vet_bound of 'Build_M_F_DD_Transiente'", 12, sizeof(double)); //
	Ub = (double*) mycalloc("Ub of 'Build_M_F_DD_Transiente'", 4, sizeof(double));
	Ue = (double*) mycalloc("Ue of 'Build_M_F_DD_Transiente'", 12, sizeof(double));
	dUe = (double*) mycalloc("dUe of 'Build_M_F_DD_Transiente'", 12, sizeof(double));
	dUb = (double*) mycalloc("dUb of 'Build_M_F_DD_Transiente'", 4, sizeof(double));

	gradUx = (double*) mycalloc("gradUx of 'Build_M_F_DD_Transiente'", 4, sizeof(double));
	gradUy = (double*) mycalloc("gradUy of 'Build_M_F_DD_Transiente'", 4, sizeof(double));
	
	dzero(neq+1, F); // F do sistema M*Da = F*
	setzeros(0, Parameters,MatrixData);
	
	for (e = 0; e < nel; e++){
		// *** Global node that composes the element
		J1 = Element[e].Vertex[0];
		J2 = Element[e].Vertex[1];
		J3 = Element[e].Vertex[2];

		// *** Nodal coordinates and differential operator
		x1 = Node[J1].x;
		x2 = Node[J2].x;
		x3 = Node[J3].x;
		
		y1 = Node[J1].y;
		y2 = Node[J2].y;
		y3 = Node[J3].y;
		
		gamma = Gamma(third*(x1+x2+x3), third*(y1+y2+y3));

		y23 = y2 - y3;
		y31 = y3 - y1;
		y12 = y1 - y2;

		x32 = x3 - x2;
		x13 = x1 - x3;
		x21 = x2 - x1;

		// *** Element Area
		twoArea = fabs(x21 * y31 - x13 * y12);
		Area = 0.5*twoArea; 
		
		// *** Fill Ue with the values ​​of the previous solution or workaround for each degree of liberty of the element node

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



		// *** The triangle centroid
		for (i = 0; i < 4; i++){
			Ub[i] = (Ue[i] + Ue[i+4] + Ue[i+8]) * third;
			dUb[i] = (dUe[i] + dUe[i+4] + dUe[i+8]) * third;
		}

		// *** Calculation of Gradient (gradu = Bu)
		for (i = 0; i < 4; i++)
		{
			gradUx[i] = (Ue[i]*y23 + Ue[i+4]*y31 + Ue[i+8]*y12) / twoArea;
			gradUy[i] = (Ue[i]*x32 + Ue[i+4]*x13 + Ue[i+8]*x21) / twoArea;
		}

		// *** Jacobian matrices Ax and Ay : primitive variables
		double rho = Ub[0];
		double v1 = Ub[1] / Ub[0];
		double v2 = Ub[2] / Ub[0];
		double E  = Ub[3] / Ub[0];
		double norma_U23 = (Ub[1] * Ub[1]) + (Ub[2] * Ub[2]);

		/*Calculation of the local Mach number*/
		/*
		double Mach, Cstab;
		Mach = sqrt( v1*v1 + v2*v2 )/sqrt( gamma*(gamma - 1.0) * ( E - 0.5 * ( v1*v1 + v2*v2 ) ) );
				
		if(Mach>Mmax)
			Mmax = Mach;
		if(Mach<Mmin)
			Mmin = Mach;
		*/
			

		// *** Ax coefficients
		double Ax[4][4];
		
		Ax[0][0] = 0.0;
		Ax[0][1] = 1.0;
		Ax[0][2] = 0.0;
		Ax[0][3] = 0.0;
		Ax[1][0] = -(v1 * v1) + (((gamma - 1) * norma_U23 * 0.5) / (rho * rho));
		Ax[1][1] = (3.0 - gamma) * v1;
		Ax[1][2] = -(gamma - 1.0) * v2;
		Ax[1][3] = gamma - 1.0;
		Ax[2][0] = - v1*v2;
		Ax[2][1] = v2;
		Ax[2][2] = v1;
		Ax[2][3] = 0.0;
		Ax[3][0] = (-gamma * v1 * E) + ((gamma - 1.0) * v1 * (norma_U23 / (rho * rho)));
		Ax[3][1] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v1 * v1)));
		Ax[3][2] = -((gamma - 1.0) * v1 * v2);
		Ax[3][3] = gamma * v1;
	
		// *** Ay coefficients
		double Ay[4][4];

		Ay[0][0] = 0.0;
		Ay[0][1] = 0.0;
		Ay[0][2] = 1.0;
		Ay[0][0] = 0.0;
		Ay[1][0] = Ax[2][0];
		Ay[1][1] = v2;
		Ay[1][2] = v1;
		Ay[1][3] = 0.0;
		Ay[2][0] = -(v2 * v2) + (((gamma - 1.0) * norma_U23 * 0.5) / (rho * rho)); 
		Ay[2][1] = -(gamma - 1.0) * v1;
		Ay[2][2] = (3.0 - gamma) * v2;
		Ay[2][3] = Ax[1][3];
		Ay[3][0] = (- gamma * v2 * E) + (((gamma - 1.0) * v2 * norma_U23) / (rho * rho));
		Ay[3][1] = Ax[3][2];
		Ay[3][2] = (gamma * E) - ((gamma - 1.0) * 0.5 * ((norma_U23 / (rho * rho)) + (2.0 * v2 * v2)));
		Ay[3][3] = gamma * v2;
		

		/* ----- Comentado, pois a norma de A0^(-1) é usada apenas no DD		

		// *** Inversa da Metrica Riemmaniana A0^(-1)
		double rho_i = Ub[3] - ((norma_U23 * 0.5) / Ub[0]); 
		
		// Variaveis de entropia do baricentro
		double V2, V3, V4;
		V2 = Ub[1] / rho_i;
		V3 = Ub[2] / rho_i;
		V4 = (- Ub[0]) / rho_i;

		// A0^(-1) coefficients
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

		*/
		double A0[4][4];		
		delta = ShockCapture(tolerance, delta_old, gradUx, gradUy, Ax, Ay, A0, dUb, y23, y31, y12, x32, x13, x21, twoArea, e, Parameters->invY, Ub); 


		//------------------------------------------------------------------------------
		//  MONTAGENS DAS MATRIZES
		//------------------------------------------------------------------------------

		//Coeficientes utilizados na matriz Khh e KBh

		// Matriz Ka = y23Ax + x32Ay
		double Ka11, Ka12, Ka13, Ka14, Ka21, Ka22, Ka23, Ka24, Ka31, Ka32, Ka33, Ka34, Ka41, Ka42, Ka43, Ka44;
		Ka11 = y23 * Ax[0][0] + x32 * Ay[0][0];
		Ka12 = y23 * Ax[0][1] + x32 * Ay[0][1];
		Ka13 = y23 * Ax[0][2] + x32 * Ay[0][2];
		Ka14 = y23 * Ax[0][3] + x32 * Ay[0][3];
		Ka21 = y23 * Ax[1][0] + x32 * Ay[1][0];
		Ka22 = y23 * Ax[1][1] + x32 * Ay[1][1];
		Ka23 = y23 * Ax[1][2] + x32 * Ay[1][2];
		Ka24 = y23 * Ax[1][3] + x32 * Ay[1][3];
		Ka31 = y23 * Ax[2][0] + x32 * Ay[2][0];
		Ka32 = y23 * Ax[2][1] + x32 * Ay[2][1];
		Ka33 = y23 * Ax[2][2] + x32 * Ay[2][2];
		Ka34 = y23 * Ax[2][3] + x32 * Ay[2][3];
		Ka41 = y23 * Ax[3][0] + x32 * Ay[3][0];
		Ka42 = y23 * Ax[3][1] + x32 * Ay[3][1];
		Ka43 = y23 * Ax[3][2] + x32 * Ay[3][2];
		Ka44 = y23 * Ax[3][3] + x32 * Ay[3][3];

		// Matriz Kb = y31Ax + x13Ay
		double Kb11, Kb12, Kb13, Kb14, Kb21, Kb22, Kb23, Kb24, Kb31, Kb32, Kb33, Kb34, Kb41, Kb42, Kb43, Kb44;
		Kb11 = y31 * Ax[0][0] + x13 * Ay[0][0];
		Kb12 = y31 * Ax[0][1] + x13 * Ay[0][1];
		Kb13 = y31 * Ax[0][2] + x13 * Ay[0][2];
		Kb14 = y31 * Ax[0][3] + x13 * Ay[0][3];
		Kb21 = y31 * Ax[1][0] + x13 * Ay[1][0];
		Kb22 = y31 * Ax[1][1] + x13 * Ay[1][1];
		Kb23 = y31 * Ax[1][2] + x13 * Ay[1][2];
		Kb24 = y31 * Ax[1][3] + x13 * Ay[1][3];
		Kb31 = y31 * Ax[2][0] + x13 * Ay[2][0];
		Kb32 = y31 * Ax[2][1] + x13 * Ay[2][1];
		Kb33 = y31 * Ax[2][2] + x13 * Ay[2][2];
		Kb34 = y31 * Ax[2][3] + x13 * Ay[2][3];
		Kb41 = y31 * Ax[3][0] + x13 * Ay[3][0];
		Kb42 = y31 * Ax[3][1] + x13 * Ay[3][1];
		Kb43 = y31 * Ax[3][2] + x13 * Ay[3][2];
		Kb44 = y31 * Ax[3][3] + x13 * Ay[3][3];

		// Matriz Kc = y12Ax + x21Ay
		double Kc11, Kc12, Kc13, Kc14, Kc21, Kc22, Kc23, Kc24, Kc31, Kc32, Kc33, Kc34, Kc41, Kc42, Kc43, Kc44;
		Kc11 = y12 * Ax[0][0] + x21 * Ay[0][0];
		Kc12 = y12 * Ax[0][1] + x21 * Ay[0][1];
		Kc13 = y12 * Ax[0][2] + x21 * Ay[0][2];
		Kc14 = y12 * Ax[0][3] + x21 * Ay[0][3];
		Kc21 = y12 * Ax[1][0] + x21 * Ay[1][0];
		Kc22 = y12 * Ax[1][1] + x21 * Ay[1][1];
		Kc23 = y12 * Ax[1][2] + x21 * Ay[1][2];
		Kc24 = y12 * Ax[1][3] + x21 * Ay[1][3];
		Kc31 = y12 * Ax[2][0] + x21 * Ay[2][0];
		Kc32 = y12 * Ax[2][1] + x21 * Ay[2][1];
		Kc33 = y12 * Ax[2][2] + x21 * Ay[2][2];
		Kc34 = y12 * Ax[2][3] + x21 * Ay[2][3];
		Kc41 = y12 * Ax[3][0] + x21 * Ay[3][0];
		Kc42 = y12 * Ax[3][1] + x21 * Ay[3][1];
		Kc43 = y12 * Ax[3][2] + x21 * Ay[3][2];
		Kc44 = y12 * Ax[3][3] + x21 * Ay[3][3];

		// *** Matriz de Massa Mhh 12x12 que é igual a matriz de massa de Galerkin    
		double Mhh1, Mhh2;
		Mhh1 = Area * sixth;
		Mhh2 = Area / 12.0;
		
		// *** Matriz de Massa MhB 12x4 e MBh 4x12, onde os coeficientes são iguais, assim calculamos apenas MhB
		double MhB;
		MhB = (3.0 * Area) / 20.0;
		
		// *** Matriz de Massa MBB 4x4
		double MBB;
		MBB = (81.0 * Area) / 280.0;
		
		// *** Matriz de Rigidez advinda de Galerkin Khhg 12x12 onde as 4 primeiras linhas repetem duas vezes, segue calculo das 4 primeiras linhas
		double Khhg11, Khhg12, Khhg13, Khhg14, Khhg21, Khhg22, Khhg23, Khhg24, Khhg31, Khhg32, Khhg33, Khhg34, Khhg41, Khhg42, Khhg43, Khhg44;
		Khhg11 = Ka11 * sixth;
		Khhg12 = Ka12 * sixth;
		Khhg13 = Ka13 * sixth;
		Khhg14 = Ka14 * sixth;
		Khhg21 = Ka21 * sixth;
		Khhg22 = Ka22 * sixth;
		Khhg23 = Ka23 * sixth;
		Khhg24 = Ka24 * sixth;
		Khhg31 = Ka31 * sixth;
		Khhg32 = Ka32 * sixth;
		Khhg33 = Ka33 * sixth;
		Khhg34 = Ka34 * sixth;
		Khhg41 = Ka41 * sixth;
		Khhg42 = Ka42 * sixth;
		Khhg43 = Ka43 * sixth;
		Khhg44 = Ka44 * sixth;
		
		double Khhg15, Khhg16, Khhg17, Khhg18, Khhg25, Khhg26, Khhg27, Khhg28, Khhg35, Khhg36, Khhg37, Khhg38, Khhg45, Khhg46, Khhg47, Khhg48;
		Khhg15 = Kb11 * sixth;
		Khhg16 = Kb12 * sixth;
		Khhg17 = Kb13 * sixth;
		Khhg18 = Kb14 * sixth;
		Khhg25 = Kb21 * sixth;
		Khhg26 = Kb22 * sixth;
		Khhg27 = Kb23 * sixth;
		Khhg28 = Kb24 * sixth;
		Khhg35 = Kb31 * sixth;
		Khhg36 = Kb32 * sixth;
		Khhg37 = Kb33 * sixth;
		Khhg38 = Kb34 * sixth;
		Khhg45 = Kb41 * sixth;
		Khhg46 = Kb42 * sixth;
		Khhg47 = Kb43 * sixth;
		Khhg48 = Kb44 * sixth;
		
		double Khhg19, Khhg110, Khhg111, Khhg112, Khhg29, Khhg210, Khhg211, Khhg212, Khhg39, Khhg310, Khhg311, Khhg312, Khhg49, Khhg410, Khhg411, Khhg412;
		Khhg19  = Kc11 * sixth;
		Khhg110 = Kc12 * sixth;
		Khhg111 = Kc13 * sixth;
		Khhg112 = Kc14 * sixth;
		Khhg29  = Kc21 * sixth;
		Khhg210 = Kc22 * sixth;
		Khhg211 = Kc23 * sixth;
		Khhg212 = Kc24 * sixth;
		Khhg39  = Kc31 * sixth;
		Khhg310 = Kc32 * sixth;
		Khhg311 = Kc33 * sixth;
		Khhg312 = Kc34 * sixth;
		Khhg49  = Kc41 * sixth;
		Khhg410 = Kc42 * sixth;
		Khhg411 = Kc43 * sixth;
		Khhg412 = Kc44 * sixth;		

		// *** Matriz de Rigidez advinda do DD Khhdd 12x12 onde a matriz é simétrica e bloco diagonal, segue calculo das 6 componentes
		double Khhdd11, Khhdd15, Khhdd19, Khhdd55, Khhdd59, Khhdd99;
		double delta4area = delta / (4.0 * Area);
		
		Khhdd11 = delta4area * (y23 * y23 + x32 * x32);
		Khhdd15 = delta4area * (y23 * y31 + x32 * x13);
		Khhdd19 = delta4area * (y23 * y12 + x32 * x21);
		Khhdd55 = delta4area * (y31 * y31 + x13 * x13);
		Khhdd59 = delta4area * (y31 * y12 + x13 * x21);
		Khhdd99 = delta4area * (y12 * y12 + x21 * x21);
		
		// *** Matriz de Rigidez KhB 12x4
		double a1 = ninefortieth * (   y23  + 2.0*y31 + 2.0*y12);
		double a2 = ninefortieth * (   x32  + 2.0*x13 + 2.0*x21);
		double b1 = ninefortieth * (2.0*y23 +   y31   + 2.0*y12);
		double b2 = ninefortieth * (2.0*x32 +   x13   + 2.0*x21);
		double c1 = ninefortieth * (2.0*y23 + 2.0*y31 +   y12  );
		double c2 = ninefortieth * (2.0*x32 + 2.0*x13 +   x21  );
		
		double KhB11, KhB12, KhB13, KhB14, KhB21, KhB22, KhB23, KhB24, KhB31, KhB32, KhB33, KhB34, KhB41, KhB42, KhB43, KhB44;
		
		KhB11 = a1*Ax[0][0] + a2*Ay[0][0];
		KhB12 = a1*Ax[0][1] + a2*Ay[0][1];
		KhB13 = a1*Ax[0][2] + a2*Ay[0][2];
		KhB14 = a1*Ax[0][3] + a2*Ay[0][3];
		KhB21 = a1*Ax[1][0] + a2*Ay[1][0];
		KhB22 = a1*Ax[1][1] + a2*Ay[1][1];		
		KhB23 = a1*Ax[1][2] + a2*Ay[1][2];
		KhB24 = a1*Ax[1][3] + a2*Ay[1][3];
		KhB31 = a1*Ax[2][0] + a2*Ay[2][0];
		KhB32 = a1*Ax[2][1] + a2*Ay[2][1];
		KhB33 = a1*Ax[2][2] + a2*Ay[2][2];
		KhB34 = a1*Ax[2][3] + a2*Ay[2][3];
		KhB41 = a1*Ax[3][0] + a2*Ay[3][0];
		KhB42 = a1*Ax[3][1] + a2*Ay[3][1];
		KhB43 = a1*Ax[3][2] + a2*Ay[3][2];
		KhB44 = a1*Ax[3][3] + a2*Ay[3][3];
		
		double KhB51, KhB52, KhB53, KhB54, KhB61, KhB62, KhB63, KhB64, KhB71, KhB72, KhB73, KhB74, KhB81, KhB82, KhB83, KhB84;
		
		KhB51 = b1*Ax[0][0] + b2*Ay[0][0];
		KhB52 = b1*Ax[0][1] + b2*Ay[0][1];
		KhB53 = b1*Ax[0][2] + b2*Ay[0][2];
		KhB54 = b1*Ax[0][3] + b2*Ay[0][3];
		KhB61 = b1*Ax[1][0] + b2*Ay[1][0];
		KhB62 = b1*Ax[1][1] + b2*Ay[1][1];		
		KhB63 = b1*Ax[1][2] + b2*Ay[1][2];
		KhB64 = b1*Ax[1][3] + b2*Ay[1][3];
		KhB71 = b1*Ax[2][0] + b2*Ay[2][0];
		KhB72 = b1*Ax[2][1] + b2*Ay[2][1];
		KhB73 = b1*Ax[2][2] + b2*Ay[2][2];
		KhB74 = b1*Ax[2][3] + b2*Ay[2][3];
		KhB81 = b1*Ax[3][0] + b2*Ay[3][0];
		KhB82 = b1*Ax[3][1] + b2*Ay[3][1];
		KhB83 = b1*Ax[3][2] + b2*Ay[3][2];
		KhB84 = b1*Ax[3][3] + b2*Ay[3][3];
		
		double KhB91, KhB92, KhB93, KhB94, KhB101, KhB102, KhB103, KhB104, KhB111, KhB112, KhB113, KhB114, KhB121, KhB122, KhB123, KhB124;
		
		KhB91 = c1*Ax[0][0] + c2*Ay[0][0];
		KhB92 = c1*Ax[0][1] + c2*Ay[0][1];
		KhB93 = c1*Ax[0][2] + c2*Ay[0][2];
		KhB94 = c1*Ax[0][3] + c2*Ay[0][3];
		KhB101 = c1*Ax[1][0] + c2*Ay[1][0];
		KhB102 = c1*Ax[1][1] + c2*Ay[1][1];		
		KhB103 = c1*Ax[1][2] + c2*Ay[1][2];
		KhB104 = c1*Ax[1][3] + c2*Ay[1][3];
		KhB111 = c1*Ax[2][0] + c2*Ay[2][0];
		KhB112 = c1*Ax[2][1] + c2*Ay[2][1];
		KhB113 = c1*Ax[2][2] + c2*Ay[2][2];
		KhB114 = c1*Ax[2][3] + c2*Ay[2][3];
		KhB121 = c1*Ax[3][0] + c2*Ay[3][0];
		KhB122 = c1*Ax[3][1] + c2*Ay[3][1];
		KhB123 = c1*Ax[3][2] + c2*Ay[3][2];
		KhB124 = c1*Ax[3][3] + c2*Ay[3][3];
		
		// *** Matriz de Rigidez KBh 4x12
		double KBh11, KBh12, KBh13, KBh14, KBh21, KBh22, KBh23, KBh24, KBh31, KBh32, KBh33, KBh34, KBh41, KBh42, KBh43, KBh44;
		KBh11 = ninefortieth * Ka11;
		KBh12 = ninefortieth * Ka12;
		KBh13 = ninefortieth * Ka13;
		KBh14 = ninefortieth * Ka14;
		KBh21 = ninefortieth * Ka21;
		KBh22 = ninefortieth * Ka22;
		KBh23 = ninefortieth * Ka23;
		KBh24 = ninefortieth * Ka24;
		KBh31 = ninefortieth * Ka31;
		KBh32 = ninefortieth * Ka32;
		KBh33 = ninefortieth * Ka33;
		KBh34 = ninefortieth * Ka34;
		KBh41 = ninefortieth * Ka41;
		KBh42 = ninefortieth * Ka42;
		KBh43 = ninefortieth * Ka43;
		KBh44 = ninefortieth * Ka44;
		
		double KBh15, KBh16, KBh17, KBh18, KBh25, KBh26, KBh27, KBh28, KBh35, KBh36, KBh37, KBh38, KBh45, KBh46, KBh47, KBh48;
		KBh15 = ninefortieth * Kb11;
		KBh16 = ninefortieth * Kb12;
		KBh17 = ninefortieth * Kb13;
		KBh18 = ninefortieth * Kb14;
		KBh25 = ninefortieth * Kb21;
		KBh26 = ninefortieth * Kb22;
		KBh27 = ninefortieth * Kb23;
		KBh28 = ninefortieth * Kb24;
		KBh35 = ninefortieth * Kb31;
		KBh36 = ninefortieth * Kb32;
		KBh37 = ninefortieth * Kb33;
		KBh38 = ninefortieth * Kb34;
		KBh45 = ninefortieth * Kb41;
		KBh46 = ninefortieth * Kb42;
		KBh47 = ninefortieth * Kb43;
		KBh48 = ninefortieth * Kb44;
		
		double KBh19, KBh110, KBh111, KBh112, KBh29, KBh210, KBh211, KBh212, KBh39, KBh310, KBh311, KBh312, KBh49, KBh410, KBh411, KBh412;
		KBh19 = ninefortieth * Kc11;
		KBh110 = ninefortieth * Kc12;
		KBh111 = ninefortieth * Kc13;
		KBh112 = ninefortieth * Kc14;
		KBh29 = ninefortieth * Kc21;
		KBh210 = ninefortieth * Kc22;
		KBh211 = ninefortieth * Kc23;
		KBh212 = ninefortieth * Kc24;
		KBh39 = ninefortieth * Kc31;
		KBh310 = ninefortieth * Kc32;
		KBh311 = ninefortieth * Kc33;
		KBh312 = ninefortieth * Kc34;
		KBh49 = ninefortieth * Kc41;
		KBh410 = ninefortieth * Kc42;
		KBh411 = ninefortieth * Kc43;
		KBh412 = ninefortieth * Kc44;
		

		delta = Delta_YZBetaNMV(tolerance, delta_old_NMV, gradUx, gradUy, Ax, Ay, A0, dUb, y23, y31, y12, x32, x13, x21, twoArea, e, Parameters->invY, Ub);
		
		// *** Matriz de Rigidez KBB 4x4
		double KBB;
		
		//delta = delta/Cstab;
		KBB = (81.0*delta*(y23*y23 + y31*y31 + y12*y12 + y23*y31 + y23*y12 + y31*y12 + x32*x32 + x13*x13 + x21*x21 + x32*x13 + x32*x21 + x13*x21))/(40.0*Area);
		
		// Calculo auxiliar
		double alpha_dt = alpha*delta_t;
		
		// M1 = Mhh + alpha*dt*Khh
		double M1[12][12];
		// no1 - no1 
		M1[0][0] = Mhh1 + alpha_dt*(Khhg11 + Khhdd11);
		M1[0][1] =      + alpha_dt*Khhg12;
		M1[0][2] =      + alpha_dt*Khhg13;
		M1[0][3] =      + alpha_dt*Khhg14;
		M1[1][0] =      + alpha_dt*Khhg21;
		M1[1][1] = Mhh1 + alpha_dt*(Khhg22 + Khhdd11);
		M1[1][2] =      + alpha_dt*Khhg23;
		M1[1][3] =      + alpha_dt*Khhg24;
		M1[2][0] =      + alpha_dt*Khhg31;
		M1[2][1] =      + alpha_dt*Khhg32;
		M1[2][2] = Mhh1 + alpha_dt*(Khhg33 + Khhdd11);
		M1[2][3] =      + alpha_dt*Khhg34;
		M1[3][0] =      + alpha_dt*Khhg41;
		M1[3][1] =      + alpha_dt*Khhg42;
		M1[3][2] =      + alpha_dt*Khhg43;
		M1[3][3] = Mhh1 + alpha_dt*(Khhg44 + Khhdd11);
		
		// no1 - no2 
		M1[0][4] = Mhh2 + alpha_dt*(Khhg15 + Khhdd15);
		M1[0][5] =      + alpha_dt*Khhg16;
		M1[0][6] =      + alpha_dt*Khhg17;
		M1[0][7] =      + alpha_dt*Khhg18;
		M1[1][4] =      + alpha_dt*Khhg25;
		M1[1][5] = Mhh2 + alpha_dt*(Khhg26 + Khhdd15);
		M1[1][6] =      + alpha_dt*Khhg27;
		M1[1][7] =      + alpha_dt*Khhg28;
		M1[2][4] =      + alpha_dt*Khhg35;
		M1[2][5] =      + alpha_dt*Khhg36;
		M1[2][6] = Mhh2 + alpha_dt*(Khhg37 + Khhdd15);
		M1[2][7] =      + alpha_dt*Khhg38;
		M1[3][4] =      + alpha_dt*Khhg45;
		M1[3][5] =      + alpha_dt*Khhg46;
		M1[3][6] =      + alpha_dt*Khhg47;
		M1[3][7] = Mhh2 + alpha_dt*(Khhg48 + Khhdd15);
		
		// no1 - no3 
		M1[0][8]  = Mhh2 + alpha_dt*(Khhg19 + Khhdd19); 
		M1[0][9]  =      + alpha_dt*Khhg110;
		M1[0][10] =      + alpha_dt*Khhg111;
		M1[0][11] =      + alpha_dt*Khhg112;
		M1[1][8]  =      + alpha_dt*Khhg29;
		M1[1][9]  = Mhh2 + alpha_dt*(Khhg210 + Khhdd19);
		M1[1][10] =      + alpha_dt*Khhg211;
		M1[1][11] =      + alpha_dt*Khhg212;
		M1[2][8]  =      + alpha_dt*Khhg39;
		M1[2][9]  =      + alpha_dt*Khhg310;
		M1[2][10] = Mhh2 + alpha_dt*(Khhg311 + Khhdd19);
		M1[2][11] =      + alpha_dt*Khhg312;
		M1[3][8]  =      + alpha_dt*Khhg49;
		M1[3][9]  =      + alpha_dt*Khhg410;
		M1[3][10] =      + alpha_dt*Khhg411;
		M1[3][11] = Mhh2 + alpha_dt*(Khhg412 + Khhdd19);
		
		// no2 - no1 
		M1[4][0] = Mhh2 + alpha_dt*(Khhg11 + Khhdd15);
		M1[4][1] =      + alpha_dt*Khhg12;
		M1[4][2] =      + alpha_dt*Khhg13;
		M1[4][3] =      + alpha_dt*Khhg14;
		M1[5][0] =      + alpha_dt*Khhg21;
		M1[5][1] = Mhh2 + alpha_dt*(Khhg22 + Khhdd15);
		M1[5][2] =      + alpha_dt*Khhg23;
		M1[5][3] =      + alpha_dt*Khhg24;
		M1[6][0] =      + alpha_dt*Khhg31;
		M1[6][1] =      + alpha_dt*Khhg32;
		M1[6][2] = Mhh2 + alpha_dt*(Khhg33 + Khhdd15);
		M1[6][3] =      + alpha_dt*Khhg34;
		M1[7][0] =      + alpha_dt*Khhg41;
		M1[7][1] =      + alpha_dt*Khhg42;
		M1[7][2] =      + alpha_dt*Khhg43;
		M1[7][3] = Mhh2 + alpha_dt*(Khhg44 + Khhdd15);
		
		// no2 - no2 
		M1[4][4] = Mhh1 + alpha_dt*(Khhg15 + Khhdd55);
		M1[4][5] =      + alpha_dt*Khhg16;
		M1[4][6] =      + alpha_dt*Khhg17;
		M1[4][7] =      + alpha_dt*Khhg18;
		M1[5][4] =      + alpha_dt*Khhg25;
		M1[5][5] = Mhh1 + alpha_dt*(Khhg26 + Khhdd55);
		M1[5][6] =      + alpha_dt*Khhg27;
		M1[5][7] =      + alpha_dt*Khhg28;
		M1[6][4] =      + alpha_dt*Khhg35;
		M1[6][5] =      + alpha_dt*Khhg36;
		M1[6][6] = Mhh1 + alpha_dt*(Khhg37 + Khhdd55);
		M1[6][7] =      + alpha_dt*Khhg38;
		M1[7][4] =      + alpha_dt*Khhg45;
		M1[7][5] =      + alpha_dt*Khhg46;
		M1[7][6] =      + alpha_dt*Khhg47;
		M1[7][7] = Mhh1 + alpha_dt*(Khhg48 + Khhdd55);
		
		// no2 - no3 
		M1[4][8]  = Mhh2 + alpha_dt*(Khhg19 + Khhdd59); 
		M1[4][9]  =      + alpha_dt*Khhg110;
		M1[4][10] =      + alpha_dt*Khhg111;
		M1[4][11] =      + alpha_dt*Khhg112;
		M1[5][8]  =      + alpha_dt*Khhg29;
		M1[5][9]  = Mhh2 + alpha_dt*(Khhg210 + Khhdd59);
		M1[5][10] =      + alpha_dt*Khhg211;
		M1[5][11] =      + alpha_dt*Khhg212;
		M1[6][8]  =      + alpha_dt*Khhg39;
		M1[6][9]  =      + alpha_dt*Khhg310;
		M1[6][10] = Mhh2 + alpha_dt*(Khhg311 + Khhdd59);
		M1[6][11] =      + alpha_dt*Khhg312;
		M1[7][8]  =      + alpha_dt*Khhg49;
		M1[7][9]  =      + alpha_dt*Khhg410;
		M1[7][10] =      + alpha_dt*Khhg411;
		M1[7][11] = Mhh2 + alpha_dt*(Khhg412 + Khhdd59);
		
		// no3 - no1 
		M1[8][0]  = Mhh2 + alpha_dt*(Khhg11 + Khhdd19);
		M1[8][1]  =      + alpha_dt*Khhg12;
		M1[8][2]  =      + alpha_dt*Khhg13;
		M1[8][3]  =      + alpha_dt*Khhg14;
		M1[9][0]  =      + alpha_dt*Khhg21;
		M1[9][1]  = Mhh2 + alpha_dt*(Khhg22 + Khhdd19);
		M1[9][2]  =      + alpha_dt*Khhg23;
		M1[9][3]  =      + alpha_dt*Khhg24;
		M1[10][0] =      + alpha_dt*Khhg31;
		M1[10][1] =      + alpha_dt*Khhg32;
		M1[10][2] = Mhh2 + alpha_dt*(Khhg33 + Khhdd19);
		M1[10][3] =      + alpha_dt*Khhg34;
		M1[11][0] =      + alpha_dt*Khhg41;
		M1[11][1] =      + alpha_dt*Khhg42;
		M1[11][2] =      + alpha_dt*Khhg43;
		M1[11][3] = Mhh2 + alpha_dt*(Khhg44 + Khhdd19);
		
		// no3 - no2 
		M1[8][4]  = Mhh2 + alpha_dt*(Khhg15 + Khhdd59);
		M1[8][5]  =      + alpha_dt*Khhg16;
		M1[8][6]  =      + alpha_dt*Khhg17;
		M1[8][7]  =      + alpha_dt*Khhg18;
		M1[9][4]  =      + alpha_dt*Khhg25;
		M1[9][5]  = Mhh2 + alpha_dt*(Khhg26 + Khhdd59);
		M1[9][6]  =      + alpha_dt*Khhg27;
		M1[9][7]  =      + alpha_dt*Khhg28;
		M1[10][4] =      + alpha_dt*Khhg35;
		M1[10][5] =      + alpha_dt*Khhg36;
		M1[10][6] = Mhh2 + alpha_dt*(Khhg37 + Khhdd59);
		M1[10][7] =      + alpha_dt*Khhg38;
		M1[11][4] =      + alpha_dt*Khhg45;
		M1[11][5] =      + alpha_dt*Khhg46;
		M1[11][6] =      + alpha_dt*Khhg47;
		M1[11][7] = Mhh2 + alpha_dt*(Khhg48 + Khhdd59);
		
		// no3 - no3 
		M1[8][8]   = Mhh1 + alpha_dt*(Khhg19 + Khhdd99); 
		M1[8][9]   =      + alpha_dt*Khhg110;
		M1[8][10]  =      + alpha_dt*Khhg111;
		M1[8][11]  =      + alpha_dt*Khhg112;
		M1[9][8]   =      + alpha_dt*Khhg29;
		M1[9][9]   = Mhh1 + alpha_dt*(Khhg210 + Khhdd99);
		M1[9][10]  =      + alpha_dt*Khhg211;
		M1[9][11]  =      + alpha_dt*Khhg212;
		M1[10][8]  =      + alpha_dt*Khhg39; 
		M1[10][9]  =      + alpha_dt*Khhg310;
		M1[10][10] = Mhh1 + alpha_dt*(Khhg311 + Khhdd99);
		M1[10][11] =      + alpha_dt*Khhg312;
		M1[11][8]  =      + alpha_dt*Khhg49;
		M1[11][9]  =      + alpha_dt*Khhg410;
		M1[11][10] =      + alpha_dt*Khhg411;
		M1[11][11] = Mhh1 + alpha_dt*(Khhg412 + Khhdd99);
		
		// N1 = MhB + alpha*dt*KhB
		double N1[12][4];
		
		N1[0][0] = MhB + alpha_dt*KhB11;
		N1[0][1] =       alpha_dt*KhB12;
		N1[0][2] =       alpha_dt*KhB13;
		N1[0][3] =       alpha_dt*KhB14;
		N1[1][0] =       alpha_dt*KhB21;
		N1[1][1] = MhB + alpha_dt*KhB22;
		N1[1][2] =       alpha_dt*KhB23;
		N1[1][3] =       alpha_dt*KhB24; 
		N1[2][0] =       alpha_dt*KhB31;
		N1[2][1] =       alpha_dt*KhB32;
		N1[2][2] = MhB + alpha_dt*KhB33;
		N1[2][3] =       alpha_dt*KhB34;
		N1[3][0] =       alpha_dt*KhB41;
		N1[3][1] =       alpha_dt*KhB42;
		N1[3][2] =       alpha_dt*KhB43;
		N1[3][3] = MhB + alpha_dt*KhB44;
		
		N1[4][0] = MhB + alpha_dt*KhB51;
		N1[4][1] =       alpha_dt*KhB52;
		N1[4][2] =       alpha_dt*KhB53;
		N1[4][3] =       alpha_dt*KhB54;
		N1[5][0] =       alpha_dt*KhB61;
		N1[5][1] = MhB + alpha_dt*KhB62;
		N1[5][2] =       alpha_dt*KhB63;
		N1[5][3] =       alpha_dt*KhB64;
		N1[6][0] =       alpha_dt*KhB71;
		N1[6][1] =       alpha_dt*KhB72;
		N1[6][2] = MhB + alpha_dt*KhB73;
		N1[6][3] =       alpha_dt*KhB74;
		N1[7][0] =       alpha_dt*KhB81;
		N1[7][1] =       alpha_dt*KhB82;
		N1[7][2] =       alpha_dt*KhB83;
		N1[7][3] = MhB + alpha_dt*KhB84;
		
		N1[8][0]  = MhB + alpha_dt*KhB91;
		N1[8][1]  =       alpha_dt*KhB92;
		N1[8][2]  =       alpha_dt*KhB93;
		N1[8][3]  =       alpha_dt*KhB94;
		N1[9][0]  =       alpha_dt*KhB101;
		N1[9][1]  = MhB + alpha_dt*KhB102;
		N1[9][2]  =       alpha_dt*KhB103;
		N1[9][3]  =       alpha_dt*KhB104;
		N1[10][0] =       alpha_dt*KhB111;
		N1[10][1] =       alpha_dt*KhB112;
		N1[10][2] = MhB + alpha_dt*KhB113;
		N1[10][3] =       alpha_dt*KhB114;
		N1[11][0] =       alpha_dt*KhB121;
		N1[11][1] =       alpha_dt*KhB122;
		N1[11][2] =       alpha_dt*KhB123;
		N1[11][3] = MhB + alpha_dt*KhB124;
		
		//M2 = MBh + alpha*dt*KBh
		double M2[4][12];
		
		M2_out[e][0]  = M2[0][0] = MhB + alpha_dt*KBh11;
		M2_out[e][1]  = M2[0][1] =       alpha_dt*KBh12;
		M2_out[e][2]  = M2[0][2] =       alpha_dt*KBh13;
		M2_out[e][3]  = M2[0][3] =       alpha_dt*KBh14;
		M2_out[e][12] = M2[1][0] =       alpha_dt*KBh21;
		M2_out[e][13] = M2[1][1] = MhB + alpha_dt*KBh22;
		M2_out[e][14] = M2[1][2] =       alpha_dt*KBh23;
		M2_out[e][15] = M2[1][3] =       alpha_dt*KBh24; 
		M2_out[e][24] = M2[2][0] =       alpha_dt*KBh31;
		M2_out[e][25] = M2[2][1] =       alpha_dt*KBh32;
		M2_out[e][26] = M2[2][2] = MhB + alpha_dt*KBh33;
		M2_out[e][27] = M2[2][3] =       alpha_dt*KBh34;
		M2_out[e][36] = M2[3][0] =       alpha_dt*KBh41;
		M2_out[e][37] = M2[3][1] =       alpha_dt*KBh42;
		M2_out[e][38] = M2[3][2] =       alpha_dt*KBh43;
		M2_out[e][39] = M2[3][3] = MhB + alpha_dt*KBh44;
		
		M2_out[e][4]  = M2[0][4] = MhB + alpha_dt*KBh15;
		M2_out[e][5]  = M2[0][5] =       alpha_dt*KBh16;
		M2_out[e][6]  = M2[0][6] =       alpha_dt*KBh17;
		M2_out[e][7]  = M2[0][7] =       alpha_dt*KBh18;
		M2_out[e][16] = M2[1][4] =       alpha_dt*KBh25;
		M2_out[e][17] = M2[1][5] = MhB + alpha_dt*KBh26;
		M2_out[e][18] = M2[1][6] =       alpha_dt*KBh27;
		M2_out[e][19] = M2[1][7] =       alpha_dt*KBh28;
		M2_out[e][28] = M2[2][4] =       alpha_dt*KBh35;
		M2_out[e][29] = M2[2][5] =       alpha_dt*KBh36;
		M2_out[e][30] = M2[2][6] = MhB + alpha_dt*KBh37;
		M2_out[e][31] = M2[2][7] =       alpha_dt*KBh38;
		M2_out[e][40] = M2[3][4] =       alpha_dt*KBh45;
		M2_out[e][41] = M2[3][5] =       alpha_dt*KBh46;
		M2_out[e][42] = M2[3][6] =       alpha_dt*KBh47;
		M2_out[e][43] = M2[3][7] = MhB + alpha_dt*KBh48;
		 
		M2_out[e][8]  = M2[0][8]  = MhB + alpha_dt*KBh19;
		M2_out[e][9]  = M2[0][9]  =       alpha_dt*KBh110;
		M2_out[e][10] = M2[0][10] =       alpha_dt*KBh111;
		M2_out[e][11] = M2[0][11] =       alpha_dt*KBh112;
		M2_out[e][20] = M2[1][8]  =       alpha_dt*KBh29;
		M2_out[e][21] = M2[1][9]  = MhB + alpha_dt*KBh210;
		M2_out[e][22] = M2[1][10] =       alpha_dt*KBh211;
		M2_out[e][23] = M2[1][11] =       alpha_dt*KBh212;
		M2_out[e][32] = M2[2][8]  =       alpha_dt*KBh39;
		M2_out[e][33] = M2[2][9]  =       alpha_dt*KBh310;
		M2_out[e][34] = M2[2][10] = MhB + alpha_dt*KBh311;
		M2_out[e][35] = M2[2][11] =       alpha_dt*KBh312;
		M2_out[e][44] = M2[3][8]  =       alpha_dt*KBh49;
		M2_out[e][45] = M2[3][9]  =       alpha_dt*KBh410;
		M2_out[e][46] = M2[3][10] =       alpha_dt*KBh411;
		M2_out[e][47] = M2[3][11] = MhB + alpha_dt*KBh412;
		
		// N2^{-1} = (MBB + alpha*dt*KBB)^{-1}
		double invN2 = 1.0 / (MBB + alpha_dt*KBB);
		invN2_out[e] = invN2;
		
		//Montagem da matriz Me
		int k, size = NDOF*NNOEL;
		double soma;
		for(i = 0; i < size; i++){
			for(j = 0; j < size; j++){
				soma = 0.0;
				for (k = 0; k < NDOF; k++){
					soma += N1[i][k]*M2[k][j];
				}
				Me[i][j] = M1[i][j] - invN2*soma;
			}
		}

		// Valor de uB e duB no elemento
		int eNDOF = e*NDOF;
		double uBaux[4], duBaux[4];
		uBaux[0] = uB[eNDOF];
		uBaux[1] = uB[eNDOF + 1];
		uBaux[2] = uB[eNDOF + 2];
		uBaux[3] = uB[eNDOF + 3];
		
		duBaux[0] = duB[eNDOF];
		duBaux[1] = duB[eNDOF + 1];
		duBaux[2] = duB[eNDOF + 2];
		duBaux[3] = duB[eNDOF + 3];
		
		//Residuo
		double R1[12], R2[4];
		
		// R1 = Fh - (Mhh*du + MhB*duB) - (Khh*u + KhB*uB)
		// Primeiro a parte da submalha: como Fh = 0, temos R1 = - MhB*duB - KhB*uB;
		R1[0]  = - MhB*duBaux[0] - (KhB11*uBaux[0]  + KhB12*uBaux[1]  + KhB13*uBaux[2]  + KhB14*uBaux[3]);
		R1[1]  = - MhB*duBaux[1] - (KhB21*uBaux[0]  + KhB22*uBaux[1]  + KhB23*uBaux[2]  + KhB24*uBaux[3]);
		R1[2]  = - MhB*duBaux[2] - (KhB31*uBaux[0]  + KhB32*uBaux[1]  + KhB33*uBaux[2]  + KhB34*uBaux[3]);
		R1[3]  = - MhB*duBaux[3] - (KhB41*uBaux[0]  + KhB42*uBaux[1]  + KhB43*uBaux[2]  + KhB44*uBaux[3]);
		R1[4]  = - MhB*duBaux[0] - (KhB51*uBaux[0]  + KhB52*uBaux[1]  + KhB53*uBaux[2]  + KhB54*uBaux[3]);
		R1[5]  = - MhB*duBaux[1] - (KhB61*uBaux[0]  + KhB62*uBaux[1]  + KhB63*uBaux[2]  + KhB64*uBaux[3]);
		R1[6]  = - MhB*duBaux[2] - (KhB71*uBaux[0]  + KhB72*uBaux[1]  + KhB73*uBaux[2]  + KhB74*uBaux[3]);
		R1[7]  = - MhB*duBaux[3] - (KhB81*uBaux[0]  + KhB82*uBaux[1]  + KhB83*uBaux[2]  + KhB84*uBaux[3]);
		R1[8]  = - MhB*duBaux[0] - (KhB91*uBaux[0]  + KhB92*uBaux[1]  + KhB93*uBaux[2]  + KhB94*uBaux[3]);
		R1[9]  = - MhB*duBaux[1] - (KhB101*uBaux[0] + KhB102*uBaux[1] + KhB103*uBaux[2] + KhB104*uBaux[3]);
		R1[10] = - MhB*duBaux[2] - (KhB111*uBaux[0] + KhB112*uBaux[1] + KhB113*uBaux[2] + KhB114*uBaux[3]);
		R1[11] = - MhB*duBaux[3] - (KhB121*uBaux[0] + KhB122*uBaux[1] + KhB123*uBaux[2] + KhB124*uBaux[3]);
		
		// Agora a parte resolvida: R1 = R1 - Mhh*du - Khh*u;
		double KhhUe[12];
		KhhUe[0]  = (Khhdd11+Khhg11)*Ue[0] +      Khhg12*Ue[1]      +      Khhg13*Ue[2]      +      Khhg14*Ue[3]      + (Khhdd15+Khhg15)*Ue[4] +      Khhg16*Ue[5]      +      Khhg17*Ue[6]      +      Khhg18*Ue[7]      + (Khhdd19+Khhg19)*Ue[8] +      Khhg110*Ue[9]      +      Khhg111*Ue[10]      +      Khhg112*Ue[11];
		KhhUe[1]  =      Khhg21*Ue[0]      + (Khhdd11+Khhg22)*Ue[1] +      Khhg23*Ue[2]      +      Khhg24*Ue[3]      +      Khhg25*Ue[4]      + (Khhdd15+Khhg26)*Ue[5] +      Khhg27*Ue[6]      +      Khhg28*Ue[7]      +      Khhg29*Ue[8]      + (Khhdd19+Khhg210)*Ue[9] +      Khhg211*Ue[10]      +      Khhg212*Ue[11];
		KhhUe[2]  =      Khhg31*Ue[0]      +      Khhg32*Ue[1]      + (Khhdd11+Khhg33)*Ue[2] +      Khhg34*Ue[3]      +      Khhg35*Ue[4]      +      Khhg36*Ue[5]      + (Khhdd15+Khhg37)*Ue[6] +      Khhg38*Ue[7]      +      Khhg39*Ue[8]      +      Khhg310*Ue[9]      + (Khhdd19+Khhg311)*Ue[10] +      Khhg312*Ue[11];
		KhhUe[3]  =      Khhg41*Ue[0]      +      Khhg42*Ue[1]      +      Khhg43*Ue[2]      + (Khhdd11+Khhg44)*Ue[3] +      Khhg45*Ue[4]      +      Khhg46*Ue[5]      +      Khhg47*Ue[6]      + (Khhdd15+Khhg48)*Ue[7] +      Khhg49*Ue[8]      +      Khhg410*Ue[9]      +      Khhg411*Ue[10]      + (Khhdd19+Khhg412)*Ue[11];
		KhhUe[4]  = (Khhdd15+Khhg11)*Ue[0] +      Khhg12*Ue[1]      +      Khhg13*Ue[2]      +      Khhg14*Ue[3]      + (Khhdd55+Khhg15)*Ue[4] +      Khhg16*Ue[5]      +      Khhg17*Ue[6]      +      Khhg18*Ue[7]      + (Khhdd59+Khhg19)*Ue[8] +      Khhg110*Ue[9]      +      Khhg111*Ue[10]      +      Khhg112*Ue[11];
		KhhUe[5]  =      Khhg21*Ue[0]      + (Khhdd15+Khhg22)*Ue[1] +      Khhg23*Ue[2]      +      Khhg24*Ue[3]      +      Khhg25*Ue[4]      + (Khhdd55+Khhg26)*Ue[5] +      Khhg27*Ue[6]      +      Khhg28*Ue[7]      +      Khhg29*Ue[8]      + (Khhdd59+Khhg210)*Ue[9] +      Khhg211*Ue[10]      +      Khhg212*Ue[11];
		KhhUe[6]  =      Khhg31*Ue[0]      +      Khhg32*Ue[1]      + (Khhdd15+Khhg33)*Ue[2] +      Khhg34*Ue[3]      +      Khhg35*Ue[4]      +      Khhg36*Ue[5]      + (Khhdd55+Khhg37)*Ue[6] +      Khhg38*Ue[7]      +      Khhg39*Ue[8]      +      Khhg310*Ue[9]      + (Khhdd59+Khhg311)*Ue[10] +      Khhg312*Ue[11];
		KhhUe[7]  =      Khhg41*Ue[0]      +      Khhg42*Ue[1]      +      Khhg43*Ue[2]      + (Khhdd15+Khhg44)*Ue[3] +      Khhg45*Ue[4]      +      Khhg46*Ue[5]      +      Khhg47*Ue[6]      + (Khhdd55+Khhg48)*Ue[7] +      Khhg49*Ue[8]      +      Khhg410*Ue[9]      +      Khhg411*Ue[10]      + (Khhdd59+Khhg412)*Ue[11];
		KhhUe[8]  = (Khhdd19+Khhg11)*Ue[0] +      Khhg12*Ue[1]      +      Khhg13*Ue[2]      +      Khhg14*Ue[3]      + (Khhdd59+Khhg15)*Ue[4] +      Khhg16*Ue[5]      +      Khhg17*Ue[6]      +      Khhg18*Ue[7]      + (Khhdd99+Khhg19)*Ue[8] +      Khhg110*Ue[9]      +      Khhg111*Ue[10]      +      Khhg112*Ue[11];
		KhhUe[9]  =      Khhg21*Ue[0]      + (Khhdd19+Khhg22)*Ue[1] +      Khhg23*Ue[2]      +      Khhg24*Ue[3]      +      Khhg25*Ue[4]      + (Khhdd59+Khhg26)*Ue[5] +      Khhg27*Ue[6]      +      Khhg28*Ue[7]      +      Khhg29*Ue[8]      + (Khhdd99+Khhg210)*Ue[9] +      Khhg211*Ue[10]      +      Khhg212*Ue[11];
		KhhUe[10] =      Khhg31*Ue[0]      +      Khhg32*Ue[1]      + (Khhdd19+Khhg33)*Ue[2] +      Khhg34*Ue[3]      +      Khhg35*Ue[4]      +      Khhg36*Ue[5]      + (Khhdd59+Khhg37)*Ue[6] +      Khhg38*Ue[7]      +      Khhg39*Ue[8]      +      Khhg310*Ue[9]      + (Khhdd99+Khhg311)*Ue[10] +      Khhg312*Ue[11]; 
		KhhUe[11] =      Khhg41*Ue[0]      +      Khhg42*Ue[1]      +      Khhg43*Ue[2]      + (Khhdd19+Khhg44)*Ue[3] +      Khhg45*Ue[4]      +      Khhg46*Ue[5]      +      Khhg47*Ue[6]      + (Khhdd59+Khhg48)*Ue[7] +      Khhg49*Ue[8]      +      Khhg410*Ue[9]      +      Khhg411*Ue[10]      + (Khhdd99+Khhg412)*Ue[11]; 
		
		R1[0]  += - (Mhh1*dUe[0] + Mhh2*dUe[4] + Mhh2*dUe[8])  - KhhUe[0];
		R1[1]  += - (Mhh1*dUe[1] + Mhh2*dUe[5] + Mhh2*dUe[9])  - KhhUe[1];
		R1[2]  += - (Mhh1*dUe[2] + Mhh2*dUe[6] + Mhh2*dUe[10]) - KhhUe[2];
		R1[3]  += - (Mhh1*dUe[3] + Mhh2*dUe[7] + Mhh2*dUe[11]) - KhhUe[3];
		R1[4]  += - (Mhh2*dUe[0] + Mhh1*dUe[4] + Mhh2*dUe[8])  - KhhUe[4];
		R1[5]  += - (Mhh2*dUe[1] + Mhh1*dUe[5] + Mhh2*dUe[9])  - KhhUe[5];
		R1[6]  += - (Mhh2*dUe[2] + Mhh1*dUe[6] + Mhh2*dUe[10]) - KhhUe[6];
		R1[7]  += - (Mhh2*dUe[3] + Mhh1*dUe[7] + Mhh2*dUe[11]) - KhhUe[7];
		R1[8]  += - (Mhh2*dUe[0] + Mhh2*dUe[4] + Mhh1*dUe[8])  - KhhUe[8];
		R1[9]  += - (Mhh2*dUe[1] + Mhh2*dUe[5] + Mhh1*dUe[9])  - KhhUe[9];
		R1[10] += - (Mhh2*dUe[2] + Mhh2*dUe[6] + Mhh1*dUe[10]) - KhhUe[10];
		R1[11] += - (Mhh2*dUe[3] + Mhh2*dUe[7] + Mhh1*dUe[11]) - KhhUe[11];
		
		R1_out[e][0] = R1[0];
		R1_out[e][1] = R1[1];
		R1_out[e][2] = R1[2];
		R1_out[e][3] = R1[3];
		R1_out[e][4] = R1[4];
		R1_out[e][5] = R1[5];
		R1_out[e][6] = R1[6];
		R1_out[e][7] = R1[7];
		R1_out[e][8] = R1[8];
		R1_out[e][9] = R1[9];
		R1_out[e][10] = R1[10];
		R1_out[e][11] = R1[11];
		
		//R2 = FB - (MBh*du + MBB*duB) - (KBh*u + KBB*uB)
		// Primeiro a parte da submalha: como FB = 0, temos R2 = - MBB*duB - KBB*uB;
		R2[0]  = - MBB*duBaux[0] - KBB*uBaux[0];
		R2[1]  = - MBB*duBaux[1] - KBB*uBaux[1];
		R2[2]  = - MBB*duBaux[2] - KBB*uBaux[2];
		R2[3]  = - MBB*duBaux[3] - KBB*uBaux[3];
		
		// Agora a parte resolvida: R2 = R2 - MBh*du - KBh*u;
		double KBhUe[4];
		KBhUe[0]  = KBh11*Ue[0] + KBh12*Ue[1] + KBh13*Ue[2] + KBh14*Ue[3] + KBh15*Ue[4] + KBh16*Ue[5] + KBh17*Ue[6] + KBh18*Ue[7] + KBh19*Ue[8] + KBh110*Ue[9] + KBh111*Ue[10] + KBh112*Ue[11];
		KBhUe[1]  = KBh21*Ue[0] + KBh22*Ue[1] + KBh23*Ue[2] + KBh24*Ue[3] + KBh25*Ue[4] + KBh26*Ue[5] + KBh27*Ue[6] + KBh28*Ue[7] + KBh29*Ue[8] + KBh210*Ue[9] + KBh211*Ue[10] + KBh212*Ue[11];
		KBhUe[2]  = KBh31*Ue[0] + KBh32*Ue[1] + KBh33*Ue[2] + KBh34*Ue[3] + KBh35*Ue[4] + KBh36*Ue[5] + KBh37*Ue[6] + KBh38*Ue[7] + KBh39*Ue[8] + KBh310*Ue[9] + KBh311*Ue[10] + KBh312*Ue[11];
		KBhUe[3]  = KBh41*Ue[0] + KBh42*Ue[1] + KBh43*Ue[2] + KBh44*Ue[3] + KBh45*Ue[4] + KBh46*Ue[5] + KBh47*Ue[6] + KBh48*Ue[7] + KBh49*Ue[8] + KBh410*Ue[9] + KBh411*Ue[10] + KBh412*Ue[11];
		
		R2[0]  += - (MhB*dUe[0] + MhB*dUe[4] + MhB*dUe[8])  - KBhUe[0];
		R2[1]  += - (MhB*dUe[1] + MhB*dUe[5] + MhB*dUe[9])  - KBhUe[1];
		R2[2]  += - (MhB*dUe[2] + MhB*dUe[6] + MhB*dUe[10]) - KBhUe[2];
		R2[3]  += - (MhB*dUe[3] + MhB*dUe[7] + MhB*dUe[11]) - KBhUe[3];
		
		R2_out[e][0] = R2[0];
		R2_out[e][1] = R2[1];
		R2_out[e][2] = R2[2];
		R2_out[e][3] = R2[3];
		
		// Montagem do vetor F 
		Fe[0]  = R1[0]  - invN2*(N1[0][0]*R2[0]  + N1[0][1]*R2[1]  + N1[0][2]*R2[2]  + N1[0][3]*R2[3]);
		Fe[1]  = R1[1]  - invN2*(N1[1][0]*R2[0]  + N1[1][1]*R2[1]  + N1[1][2]*R2[2]  + N1[1][3]*R2[3]);
		Fe[2]  = R1[2]  - invN2*(N1[2][0]*R2[0]  + N1[2][1]*R2[1]  + N1[2][2]*R2[2]  + N1[2][3]*R2[3]);
		Fe[3]  = R1[3]  - invN2*(N1[3][0]*R2[0]  + N1[3][1]*R2[1]  + N1[3][2]*R2[2]  + N1[3][3]*R2[3]);
		Fe[4]  = R1[4]  - invN2*(N1[4][0]*R2[0]  + N1[4][1]*R2[1]  + N1[4][2]*R2[2]  + N1[4][3]*R2[3]);
		Fe[5]  = R1[5]  - invN2*(N1[5][0]*R2[0]  + N1[5][1]*R2[1]  + N1[5][2]*R2[2]  + N1[5][3]*R2[3]);
		Fe[6]  = R1[6]  - invN2*(N1[6][0]*R2[0]  + N1[6][1]*R2[1]  + N1[6][2]*R2[2]  + N1[6][3]*R2[3]);
		Fe[7]  = R1[7]  - invN2*(N1[7][0]*R2[0]  + N1[7][1]*R2[1]  + N1[7][2]*R2[2]  + N1[7][3]*R2[3]);
		Fe[8]  = R1[8]  - invN2*(N1[8][0]*R2[0]  + N1[8][1]*R2[1]  + N1[8][2]*R2[2]  + N1[8][3]*R2[3]);
		Fe[9]  = R1[9]  - invN2*(N1[9][0]*R2[0]  + N1[9][1]*R2[1]  + N1[9][2]*R2[2]  + N1[9][3]*R2[3]);
		Fe[10] = R1[10] - invN2*(N1[10][0]*R2[0] + N1[10][1]*R2[1] + N1[10][2]*R2[2] + N1[10][3]*R2[3]);
		Fe[11] = R1[11] - invN2*(N1[11][0]*R2[0] + N1[11][1]*R2[1] + N1[11][2]*R2[2] + N1[11][3]*R2[3]);
		
		for (i = 0; i<NDOF; i++){
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

		for (i = 0; i < 12; i++){
			F[lm_aux[i]] += Fe[i];            
		}
		
		assembly(0, e, Me, MatrixData);
		
	}//for elemento
	
	//printf("\n\nMach min = %lf \t Mach max = %lf\n\n",Mmin,Mmax);
	
	// Liberacao dos espacos alocados
	free(dUe);
	free(dUb);
	free(Ue);
	free(gradUx);
	free(gradUy);
	free(Ub);
	free(vet_bound);
	
	return 0;
}// end build
