MPI_VectorUpdate(RANK, NPROC, Parameters, FemStructs, FemStructs->f);
		char filename[200];
		FILE *Out;

		sprintf(filename,"RANK_%d.dat",RANK);
		Out = fopen(filename,"w");
	if (RANK != 0){
		fprintf(Out,"AA%d_bef=[\n",RANK);
		for (i=0; i<Parameters->nnzero_bef; i++)
			fprintf(Out,"%.14lf ", MatrixData->AA_bef[i]);
		fprintf(Out,"];\n");
		fprintf(Out,"JA%d_bef=[\n",RANK);
		for (i=0; i<Parameters->nnzero_bef; i++)
			fprintf(Out,"%d ", MatrixData->JA_bef[i]);
		fprintf(Out,"];\n");
		fprintf(Out,"IA%d_bef=[\n",RANK);
		for (i=0; i<Parameters->neqaux_bef+1; i++)
			fprintf(Out,"%d ", MatrixData->IA_bef[i]);
		fprintf(Out,"];\n");
		fprintf(Out,"aux%d_bef=[\n",RANK);
		for (i=0; i<Parameters->neqaux_bef; i++)
			fprintf(Out,"%d ", MatrixData->aux_bef[i]);
		fprintf(Out,"];\n");
 
		for (i=0; i<Parameters->neqaux_bef; i++){
			fprintf(Out,"F%d_bef[%d]=%lf;\n", RANK, MatrixData->JA_bef[i], FemStructs->f[MatrixData->JA_bef[i]]);
		}
	}
	if (RANK != 2){
		fprintf(Out,"AA%d_aft=[\n",RANK);
		for (i=0; i<Parameters->nnzero_aft; i++)
			fprintf(Out,"%.14lf ", MatrixData->AA_aft[i]);
		fprintf(Out,"];\n");
		fprintf(Out,"JA%d_aft=[\n",RANK);
		for (i=0; i<Parameters->nnzero_aft; i++)
			fprintf(Out,"%d ", MatrixData->JA_aft[i]);
		fprintf(Out,"];\n");
		fprintf(Out,"IA%d_aft=[\n",RANK);
		for (i=0; i<Parameters->neqaux_aft+1; i++)
			fprintf(Out,"%d ", MatrixData->IA_aft[i]);
		fprintf(Out,"];\n");
		fprintf(Out,"aux%d_aft=[\n",RANK);
		for (i=0; i<Parameters->neqaux_aft; i++)
			fprintf(Out,"%d ", MatrixData->aux_aft[i]);
		fprintf(Out,"];\n");

		for (i=0; i<Parameters->neqaux_aft; i++){
			fprintf(Out,"F%d_aft[%d]=%lf;\n", RANK, MatrixData->JA_aft[i], FemStructs->f[MatrixData->JA_aft[i]+Parameters->neq_bef+Parameters->neq]);
		}
	}

	
	fprintf(Out,"AA%d=[\n",RANK);
	for (i=0; i<Parameters->nnzero; i++)
		fprintf(Out,"%.14lf ", MatrixData->AA[i]);
	fprintf(Out,"];\n");
	fprintf(Out,"JA%d=[\n",RANK);
	for (i=0; i<Parameters->nnzero; i++)
		fprintf(Out,"%d ", MatrixData->JA[i]);
	fprintf(Out,"];\n");
	fprintf(Out,"IA%d=[\n",RANK);
	for (i=0; i<Parameters->neq+1; i++)
		fprintf(Out,"%d ", MatrixData->IA[i]);
	fprintf(Out,"];\n");
 	
	for (i=0; i<Parameters->neq; i++){
		fprintf(Out,"F%d[%d]=%lf;\n", RANK, MatrixData->JA[i], FemStructs->f[MatrixData->JA[i]+Parameters->neq]);
	}

	for (i=0;i<Parameters->nnodes;i++){
		if (FemStructs->Node[i].Type!=0){
			fprintf(Out,"F[%d]=%lf (RANK %d)\n",FemStructs->Node[i].invP_id + 8*RANK, FemStructs->F[FemStructs->Node[i].invP_id+Parameters->neq_bef],RANK);
		}
	}

	fclose(Out);


