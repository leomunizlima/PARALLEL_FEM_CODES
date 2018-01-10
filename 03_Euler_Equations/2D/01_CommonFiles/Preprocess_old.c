#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

int Preprocess(int narg, char **arguments, MatrixDataType **MatrixData_out, double **F_out, double **u_out, int ***lm_out, int **lmaux_out, 
				NodeType **Node_out, ElementType **Element_out, ParametersType **Parameters_out, int **eqrho_out)
{
	int neq, neqrho, nnodes, nedge, nnzero, nel, I, J;
	int **lm, *lmaux, *eqrho;
	int size = NDOF*NNOEL;
	int size2 = size*size;
	double *F, *u;
	char FileName[300], label[300];
	FILE *InFile;
	NodeType *Node;
	ElementType *Element;
	ParametersType *Parameters;
	MatrixDataType *MatrixData;
	

	/* **************************************************************************************************************************** */
	//												Testing initial parameters
	/* **************************************************************************************************************************** */
	if (narg!=2)
	{
		printf("Use ./EulerEquations2D <Parameters file according README>\n");
		exit(1);
	}
	/* **************************************************************************************************************************** */


	/* **************************************************************************************************************************** */
	//					Reading parameters from problem setting file
	/* **************************************************************************************************************************** */
	Parameters = mycalloc("Parameters of 'Preprocess'", 1, sizeof(ParametersType));

	InFile = myfopen(arguments[1], "r");
	fscanf(InFile, "%s\t:%[^\n]", Parameters->Experiments, label);
	fscanf(InFile, "%s\t:%[^\n]", Parameters->ProblemTitle, label);
	fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->SolverTolerance), label);
	fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->CorrectionTolerance), label);
	fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->TimeIntegrationTolerance), label);
	fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->CoefficientTolerance), label);
	fscanf(InFile, "%d\t:%[^\n]", &(Parameters->itermax), label);
	fscanf(InFile, "%d\t:%[^\n]", &(Parameters->KrylovBasisVectorsQuantity), label);
	fscanf(InFile, "%s\t:%[^\n]", Parameters->Solver, label);
	fscanf(InFile, "%s\t:%[^\n]", Parameters->TimeIntegration, label);
	fscanf(InFile, "%s\t:%[^\n]", Parameters->Precondicionators, label);
	fscanf(InFile, "%lf\t:%[^\n]", &Parameters->Alpha, label);
	fscanf(InFile, "%lf\t:%[^\n]", &Parameters->DeltaT, label);
	fscanf(InFile, "%lf\t:%[^\n]", &Parameters->FinalT, label);
	fscanf(InFile, "%d\t:%[^\n]", &Parameters->NumberCorrection, label);
	fscanf(InFile, "%s\t:%[^\n]", Parameters->StopMulticorrection, label);
	fscanf(InFile, "%s\t:%[^\n]", Parameters->MatrixVectorProductScheme, label);
	fscanf(InFile, "%s\t:%[^\n]", Parameters->StabilizationForm, label);
	fscanf(InFile, "%s\t:%[^\n]", Parameters->ShockCapture, label);
	fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->invY[0]), label);
	fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->invY[1]), label);
	fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->invY[2]), label);
	fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->invY[3]), label);
	fscanf(InFile, "%d\t:%[^\n]", &(Parameters->nnodes), label);
	fscanf(InFile, "%d\t:%[^\n]", &(Parameters->nel), label);
	fclose(InFile);
	
	/* **************************************************************************************************************************** */


	/* **************************************************************************************************************************** */
	//						Reading nodes
	/* **************************************************************************************************************************** */
	sprintf(FileName,"../02_mesh/%s_%d_%d.dat", Parameters->ProblemTitle, Parameters->nnodes, Parameters->nel);
	InFile = myfopen(FileName, "r");
	fscanf(InFile, "%d", &nnodes);
	Node = (NodeType*) mycalloc("Node of 'Preprocess'", nnodes, sizeof(NodeType));
	for (I = 0; I < nnodes; I++)
	{
		fscanf(InFile, "%lf%lf%d%d%d%d", &(Node[I].x), &(Node[I].y), &(Node[I].rhoType), &(Node[I].v1Type), &(Node[I].v2Type), &(Node[I].eType));
	}
	Fill_ID(&neq, &neqrho, Node, nnodes);
	
	/* **************************************************************************************************************************** */


	/* **************************************************************************************************************************** */
	//           				Reading connection mesh
	/* **************************************************************************************************************************** */
	fscanf(InFile, "%d", &nel);
	Element = (ElementType*) mycalloc("Element of 'Preprocess'", nel, sizeof(ElementType));
	for (I = 0; I < nel; I++)
		fscanf(InFile, "%d%d%d%d", &(Element[I].Vertex[0]), &(Element[I].Vertex[1]), &(Element[I].Vertex[2]), &(Element[I].Type));
	fclose(InFile);
	
	/* **************************************************************************************************************************** */


	/* **************************************************************************************************************************** */
	//			          Memory allocations and Store strategies 
	/* **************************************************************************************************************************** */

	// Some variable inicializations
	nedge = 0;
	nnzero = 0;

	MatrixData = (MatrixDataType *) mycalloc("MatrixData of 'Preprocess'", 1, sizeof(MatrixDataType));
	F = (double*) mycalloc("F of 'Preprocess'", neq+1, sizeof(double));
	u = (double*) mycalloc("u of 'Preprocess'", neq+1, sizeof(double));
	lm = (int**) mycalloc("lm of 'Preprocess'", nel, sizeof(int*));
	lmaux = (int*) mycalloc("lmaux of 'Preprocess'", nel*size, sizeof(int));
	for (I = 0; I < nel; I++)
		lm[I] = &lmaux[I*size];
		
	
	// Set vector eqrho
	eqrho = (int*) mycalloc("eqrho of 'Preprocess'", neqrho, sizeof(int));
	J = 0;
	for (I = 0; I < nnodes; I++){
		if(Node[I].id[0] >= 0){ // a densidade é incógnita naquele nó
			eqrho[J] = Node[I].id[0];
			J++;
		}
	}
			

	//Configuring equation according to variables and boundary conditions
	Fill_LM(neq, nel, lm, Node, Element);
	
	
	if (strcasecmp(Parameters->TimeIntegration,"PREDICTOR1") == 0){	
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE") == 0){
		
			double **M, *Maux, **K, *Kaux;

			M = (double**) mycalloc("M of 'Preprocess'", nel, sizeof(double*));
			Maux = (double*) mycalloc("Maux of 'Preprocess'", nel*size2,sizeof(double));
			K = (double**) mycalloc("K of 'Preprocess'", nel, sizeof(double*));
			Kaux = (double*) mycalloc("Kaux of 'Preprocess'", nel*size2,sizeof(double));

			for (I = 0; I < nel; I++){
				M[I] = &Maux[I*size2];
				K[I] = &Kaux[I*size2];
			}

			MatrixData->A[0] = M;
			MatrixData->Aaux[0] = Maux;
			MatrixData->A[1] = K;
			MatrixData->Aaux[1] = Kaux;
			
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE") == 0){
			
			int **order, **EDGE_by_Element;
			double **M, **K, *Maux, *Kaux;		
	
			size2 = NDOF*NDOF;
			order = mycalloc("order of 'Preprocess'", nel, sizeof(int*));
			for (I = 0; I < nel; I++)
				order[I] = (int*) mycalloc("order of 'Preprocess'", NNOEL, sizeof(int));

			ede_Initialization(nnodes, neq, nel, Node, Element, &nedge, order, &lm, &lmaux, &EDGE_by_Element);
			M = (double**) mycalloc("M of 'Preprocess'", nedge, sizeof(double*));
			K = (double**) mycalloc("K of 'Preprocess'", nedge, sizeof(double*));
			Maux = (double*) mycalloc("Maux of 'Preprocess'", (nedge)*size2*(NNOEL+1),sizeof(double));
			Kaux = (double*) mycalloc("Kaux of 'Preprocess'", (nedge)*size2*(NNOEL+1),sizeof(double));

			for (I = 0; I < nedge; I++){
				M[I] = &Maux[I*size2*(NNOEL+1)];
				K[I] = &Kaux[I*size2*(NNOEL+1)];
			}

			MatrixData->A[0] = M;
			MatrixData->A[1] = K;
			MatrixData->Aaux[0] = Maux;
			MatrixData->Aaux[1] = Kaux;
			MatrixData->Scheme_by_Element = EDGE_by_Element;
			MatrixData->order = order;

		}
		else{
			printf("Matrix vector produto scheme not defined!\n\n");
			exit(1);
		}

	}else if (strcasecmp(Parameters->TimeIntegration,"PREDICTOR2") == 0){
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE") == 0){
		
			double **M, *Maux;

			M = (double**) mycalloc("M of 'Preprocess'", nel, sizeof(double*));
			Maux = (double*) mycalloc("Maux of 'Preprocess'", nel*size2,sizeof(double));
			
			for (I = 0; I < nel; I++){
				M[I] = &Maux[I*size2];
			}

			MatrixData->A[0] = M;
			MatrixData->Aaux[0] = Maux;
			
		}else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE") == 0){
			
			int **order, **EDGE_by_Element;
			double **M, *Maux;		
	
			size2 = NDOF*NDOF;
			order = mycalloc("order of 'Preprocess'", nel, sizeof(int*));
			for (I = 0; I < nel; I++)
				order[I] = (int*) mycalloc("order of 'Preprocess'", NNOEL, sizeof(int));

			ede_Initialization(nnodes, neq, nel, Node, Element, &nedge, order, &lm, &lmaux, &EDGE_by_Element);
			M = (double**) mycalloc("M of 'Preprocess'", nedge, sizeof(double*));
			Maux = (double*) mycalloc("Maux of 'Preprocess'", (nedge)*size2*(NNOEL+1),sizeof(double));

			for (I = 0; I < nedge; I++){
				M[I] = &Maux[I*size2*(NNOEL+1)];
			}

			MatrixData->A[0] = M;
			MatrixData->Aaux[0] = Maux;
			MatrixData->Scheme_by_Element = EDGE_by_Element;
			MatrixData->order = order;

		}
		else{
			printf("Matrix vector produto scheme not defined!\n\n");
			exit(1);
		}
	}else{
		printf("Time integration scheme not defined!\n\n");
		exit(1);
	}
	
	/* **************************************************************************************************************************** */

	Parameters->neq = neq;
	Parameters->neqrho = neqrho;
	Parameters->nedge = nedge;
	Parameters->nnzero = nnzero;
	*Node_out = Node;
	*Element_out = Element;
	*lm_out = lm;
	*lmaux_out = lmaux;
	*F_out = F;
	*u_out = u;
	*Parameters_out = Parameters;
	*MatrixData_out = MatrixData;
	*eqrho_out = eqrho;

	return 0;
}


