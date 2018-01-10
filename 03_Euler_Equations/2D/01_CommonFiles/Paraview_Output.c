#include "EulerEquations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Paraview_Output(int RANK, int NPROC, ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, nnodes, nel;
	char FileName[300];
	FILE *OutFile;
	double *rho, *v1, *v2, *e, *temp, *pres, aux, cv, gamma, X, Y;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	double *U = FemStructs->u;
	nnodes = Parameters->nnodes;
	nel = Parameters->nel;

	rho = (double*) mycalloc("rho of 'Paraview_Output'", nnodes, sizeof(double));
	v1 = (double*) mycalloc("v1 of 'Paraview_Output'", nnodes, sizeof(double));
	v2 = (double*) mycalloc("v2 of 'Paraview_Output'", nnodes, sizeof(double));
	e = (double*) mycalloc("e of 'Paraview_Output'", nnodes, sizeof(double));
	temp = (double*) mycalloc("temp of 'Paraview_Output'", nnodes, sizeof(double));
	pres = (double*) mycalloc("pres of 'Paraview_Output'", nnodes, sizeof(double));

	for (I = 0; I < nnodes; I++){
		X = Node[I].x;
		Y = Node[I].y;
		cv = FemFunctions->cv(X, Y);
		gamma = FemFunctions->gamma(X, Y);

		if (Node[I].Type[0] != 0)
			rho[I] = U[Node[I].invP_id[0]];
		else
			rho[I] = FemFunctions->rhopresc(X, Y);

		if (Node[I].Type[1] != 0)
			v1[I] = U[Node[I].invP_id[1]];
		else
			v1[I] = FemFunctions->v1presc(X, Y);

		if (Node[I].Type[2] != 0)
	   		v2[I] = U[Node[I].invP_id[2]];
	   	else
			v2[I] = FemFunctions->v2presc(X, Y);

		if (Node[I].Type[3] != 0)
			e[I] = U[Node[I].invP_id[3]];
		else
			e[I] = FemFunctions->epresc(X, Y);

		aux = e[I] - (v1[I]*v1[I] + v2[I]*v2[I])/(2.0*rho[I]);
		temp[I] = aux/(rho[I]*cv);
		pres[I] = (gamma - 1)*aux;
	}

	
	sprintf(FileName,"../03_output/%02d/RANK_%02d_%s_%s_%s_%s_%s_%s_%s_N%d_E%d.vtu", NPROC, RANK, Parameters->Experiments, Parameters->ProblemTitle, 
		Parameters->StabilizationForm, Parameters->ShockCapture,Parameters->TimeIntegration,Parameters->MatrixVectorProductScheme, 
		Parameters->Preconditioner, Parameters->nNodes, Parameters->nEl);
	OutFile = myfopen(FileName,"w");

	fprintf(OutFile,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	fprintf(OutFile,"\t<UnstructuredGrid>\n");
	fprintf(OutFile,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nnodes, nel);
	fprintf(OutFile,"\t\t\t<PointData Scalars=\"scalars\">\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", rho[I]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Temperature\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", temp[I]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
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

	free(rho);
	free(v1);
	free(v2);
	free(e);
	free(temp);
	free(pres);

	return 0;
}
