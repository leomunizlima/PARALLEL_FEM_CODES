int ebemv(ParametersType *Parameters, MatrixDataType *MatrixData, int structNum, double *P, double *Q, int **lm)
{

	int I, J, K, i1, i2, i3, j1, j2, j3;
	double p1, p2, p3, q1, q2, q3;
	double **A;
	int neq, nel, nnodes;
	int K1, K2, K3, K4, K5, K6, K7, K8, K9;
	int const1, const2, const3;

	neq = Parameters->neq;
	nel = Parameters->nel;
	nnodes = Parameters->nnodes;

	A = MatrixData->A[0];

	for (I = 0; I < neq; I++)
		Q[I] = 0;

	P[neq] = 0;
	Q[neq] = 0;

	const1 = 2*NDOF;
	const2 = 3*NDOF*NDOF;
	const3 = 6*NDOF*NDOF;

	for(I = 0; I < nel; I++){

		for(J = 0; J < NDOF; J++){
			i1 = lm[I][J];
			i2 = lm[I][J+NDOF];
			i3 = lm[I][J+2*NDOF];
		
			q1 = 0.0;
			q2 = 0.0;
			q3 = 0.0;			
			
			for(K = 0; K < NDOF; K++){
				j1 = lm[I][K];
				j2 = lm[I][K+NDOF];
				j3 = lm[I][K+2*NDOF];
			
				p1 = P[j1];
				p2 = P[j2];
				p3 = P[j3];	
				
				K1 = J*NDOF*NNOEL + K;
				K2 = K1 + NDOF;
				K3 = K1 + const1;
				K4 = K1 + const2;
				K5 = K4 + NDOF;
				K6 = K4 + const1;
				K7 = K1 + const3;
				K8 = K7 + NDOF;
				K9 = K7 + const1;
				q1 += A[I][K1]*p1 + A[I][K2]*p2 + A[I][K3]*p3;
				q2 += A[I][K4]*p1 + A[I][K5]*p2 + A[I][K6]*p3;
				q3 += A[I][K7]*p1 + A[I][K8]*p2 + A[I][K9]*p3;
				
			}

			Q[i1] += q1;
			Q[i2] += q2;
			Q[i3] += q3;

			Q[neq] = 0.0;

		}

	}

	return 0;
}
