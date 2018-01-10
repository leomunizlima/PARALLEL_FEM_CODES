#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Paraview_Output(int RANK, int NPROC, ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, t1, t2, t3, nnodes, nel;
	char FileName[300];
	FILE *OutFile;
	double *v1, *v2, *pres, X, Y, normu, normv, normp;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	double *U = FemStructs->u;
	nnodes = Parameters->nnodes;
	nel = Parameters->nel;
	
	v1 = (double*) mycalloc("v1 of 'Paraview_Output'", nnodes, sizeof(double));
	v2 = (double*) mycalloc("v2 of 'Paraview_Output'", nnodes, sizeof(double));
	pres = (double*) mycalloc("pres of 'Paraview_Output'", nnodes, sizeof(double));
	
	for (I = 0; I < nnodes; I++){
		t1 = Node[I].Type[0];
		t2 = Node[I].Type[1];
		t3 = Node[I].Type[2];
		//eq4 = Node[I].id[3];
		X = Node[I].x;	
		Y = Node[I].y;	
		//cv = CV(X, Y);
		//gamma = Gamma(X, Y);

		if (t1 > 0)
			v1[I] = U[Node[I].invP_id[0]];
		else
			v1[I] = FemFunctions->v1presc(X, Y);

		if (t2 > 0)
			v2[I] = U[Node[I].invP_id[1]];
		else
			v2[I] = FemFunctions->v2presc(X, Y);
	   	
	   	if (t3 > 0)
			pres[I] = U[Node[I].invP_id[2]];
	   	else
			pres[I] = FemFunctions->ppresc(X, Y);
	}
	
	// Calculo da norma de v e p
	normu = sqrt(ddot(nnodes, v1, v1));
	normv = sqrt(ddot(nnodes, v2, v2));
	normp = sqrt(ddot(nnodes, pres, pres));
	if (RANK==0)
		printf(" \n Normas: |v_x| = %.14lf, |v_y| = %.14lf, |p| = %.14lf\n", normu, normv, normp);
	
	// Calculo erro entre a solucao aprox e exata
	/*ev1 = errov1[0];
	ev2 = errov2[0];
	epres = erropres[0];
	for (I = 1; I < nnodes; I++){
		if(ev1<errov1[I])
			ev1 = errov1[I];	
		if(ev2<errov2[I])
			ev2 = errov2[I];
		if(epres<erropres[I]){
			epres = erropres[I];
			//printf("\n No do erro para a pressao: %d, %f, %f", I, Node[I].x, Node[I].y);		
		}	
	}
	printf(" \n Normas do erro: |Ev_x| = %.4lf, |Ev_y| = %.4lf, |Ep| = %.4lf (in Paraview_Output.c)", ev1, ev2, epres); */


/*	sprintf(FileName,"../03_output/%s_%s_%s_%s_%s_%s_N%d_E%d.vtu", Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->ShockCapture, */
/*			Parameters->TimeIntegration,Parameters->MatrixVectorProductScheme,Parameters->nnodes, Parameters->nel);*/
	sprintf(FileName,"../03_output/%02d/RANK_%02d_%s_%s_%s_%s_%s_N%d_E%d.vtu", NPROC, RANK, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, 
	Parameters->MatrixVectorProductScheme, Parameters->Preconditioner,Parameters->nNodes, Parameters->nEl);
	OutFile = myfopen(FileName,"w");

	fprintf(OutFile,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	fprintf(OutFile,"\t<UnstructuredGrid>\n");
	fprintf(OutFile,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nnodes, nel);
	fprintf(OutFile,"\t\t\t<PointData Scalars=\"scalars\" Vectors = \"Velocity\">\n");
//	fprintf(OutFile,"\t\t\t<PointData Vectors = \"Velocity\">\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf   %.12lf   %.12lf\n", v1[I], v2[I], 0.0);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
/*	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n");*/

/*	for (I = 0; I < nnodes; I++)*/
/*		fprintf(OutFile,"\t\t\t\t   %.12lf\n", rho[I]);*/

/*	fprintf(OutFile,"\t\t\t\t</DataArray>\n");*/
/*	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Temperature\" Format=\"ascii\">\n");*/

/*	for (I = 0; I < nnodes; I++)*/
/*		fprintf(OutFile,"\t\t\t\t   %.12lf\n", temp[I]);*/

/*	fprintf(OutFile,"\t\t\t\t</DataArray>\n");*/
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Pressure\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", pres[I]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n"); 
	fprintf(OutFile,"\t\t\t</PointData>\n");
	fprintf(OutFile,"\t\t\t<Points>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\t%.12lf\t%.12lf\n", Node[I].x,Node[I].y, 0.0);
		
	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Points>\n");
	fprintf(OutFile,"\t\t\t<Cells>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
		fprintf(OutFile,"\t\t\t\t   %d\t%d\t%d\n", Element[I].Vertex[0], Element[I].Vertex[1], Element[I].Vertex[2]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
		fprintf(OutFile,"\t\t\t\t   %d\n",(I+1)*3);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
			fprintf(OutFile,"\t\t\t\t   %d\n",5);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Cells>\n");
	fprintf(OutFile,"\t\t</Piece>\n");
	fprintf(OutFile,"\t</UnstructuredGrid>\n");
	fprintf(OutFile,"</VTKFile>\n");

	fclose(OutFile);

	free(v1);
	free(v2);
	free(pres);
		
	return 0;
}


