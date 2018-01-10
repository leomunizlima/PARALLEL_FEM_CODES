#include "solvers.h"
#include "preconditioners.h"

int pcg(ParametersType *Parameters, 
				MatrixDataType *MatrixData, 
				double *B, 
				double *X, 
				int **lm, 
				int (*precond)(ParametersType *, MatrixDataType *, double *, double *),
	 			int (*mv)(ParametersType *, MatrixDataType *, int , double *, double *, int **))
{
	int iter, itermax, neq;
	double *R, *P, *Q, *Z, gammaOld, gammaNew, tau, alpha, beta, r0, Tol;


	neq = Parameters->neq;

	// Some CG inicializations
	itermax = Parameters->LinearMaxIter;
	Tol = Parameters->SolverTolerance;

	/*************************************************************/
	//		Memory allocations for CG 
	/************************************************************/
	P = mycalloc("P of 'cg'", (neq+1),sizeof(double));
	Q = mycalloc("Q of 'cg'", (neq+1),sizeof(double));
	R = mycalloc("R of 'cg'", neq,sizeof(double));
	Z = mycalloc("Z of 'cg'", neq,sizeof(double));
	/************************************************************/

	dzero(neq,X);

	memcpy(R, B, neq*sizeof(double)); //R <-- B

	precond(Parameters, MatrixData, R, Z); // Preconditioning Z = invM * R

	memcpy(P, Z, neq*sizeof(double)); //P <-- Z

	gammaNew = ddot(neq, R, Z);
	
	r0 = gammaNew;

	#ifdef debug
		printf("gamma = %lf\n", gammaNew);
	#endif

	if (gammaNew < Tol*Tol*r0)
		return 0;
	
	for (iter=0; iter<itermax; iter++){
	
		mv(Parameters, MatrixData, 0, P, Q, lm); // Q <-- A*P
	
		tau = ddot(neq, P, Q);
	
		alpha = gammaNew/tau;

		daxpy(neq, alpha, P, X); // X <-- X + alpha*P

		daxpy(neq, -alpha, Q, R); // R <-- R - alpha*Q

		precond(Parameters, MatrixData, R, Z); // Preconditioning Z = invM * R

		gammaOld = gammaNew;

		gammaNew = ddot(neq, R, Z);

		#ifdef debug			
			printf("Gamma = %.15lf (iter = %d)\n", gammaNew,iter);
		#endif

		if (gammaNew < Tol*Tol*r0)	
			break;
		
		beta = gammaNew/gammaOld;

		dscal(neq, beta, P); // P <-- beta*P

		daxpy(neq, 1.0, Z, P); // P <-- P + Z	
	}

	#ifdef debug
		printf("cg iterations = %d\n",iter);
	#endif

	Parameters->iterations = iter;

	free(P);
	free(Q);
	free(R);
	free(Z);
	return 0;
}


