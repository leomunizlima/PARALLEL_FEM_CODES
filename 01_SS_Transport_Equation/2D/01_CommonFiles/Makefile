CC=mpicc
#CC=tau_cc.sh
IDIR= 
CFLAGS= -DSSTranspEquation2D -Ddebug -c -lm -lrt -Wall -Ofast -march=native 
GCFLAGS= -DSSTransEquation2D -Ddebug -lm -lrt -Wall -Ofast -march=native \
../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-I/opt/openmpi/include -L/opt/openmpi/lib -lmpi -libverbs -lrt -lnsl -lutil -lm -ldl -lm -Wl,--export-dynamic \
-lrt -lnsl -lutil -lm -ldl -L/usr/local/lib -lm -lgfortran -Wall -L ../../../00_CommonFiles/Solvers_and_Preconditioners/pardiso -lpardiso500-GNU481-X86-64 \
-lblas -llapack -lquadmath -fopenmp
LDFLAGS=
SOURCES= main.c \
Preprocess.c \
Process.c \
Postprocess.c \
Paraview_Output.c \
setProblem.c \
setMatrixVectorProductType.c \
setSolver.c \
setScaling.c \
setPreconditioner.c \
setStabilizationForm.c \
setzeros.c \
csr_Initialization.c \
ebe_assembly.c \
csr_assembly.c \
F_assembly.c \
Build_K_F_SUPG.c \
Build_K_F_DD.c \
CAU_ShockCapture.c \
CAU_DD_ShockCapture.c \
YZBeta_ShockCapture.c \
YZBeta_ShockCapture2.c \
Fill_LM.c \
Permutation_of_LM.c \
h_shock_2sqrtArea.c \
h_shock_Option1.c \
h_shock_Option2.c \
MPI_VectorUpdate.c \
COMPARE_array_in_Preprocess.c \
../BIGPUDIM/BIGPUDIM_Condutivity.c \
../BIGPUDIM/BIGPUDIM_Font.c \
../BIGPUDIM/BIGPUDIM_Reaction.c \
../BIGPUDIM/BIGPUDIM_Velocity.c \
../BIGPUDIM/BIGPUDIM_upresc.c \
../CARTOLA/CARTOLA_Condutivity.c \
../CARTOLA/CARTOLA_Font.c \
../CARTOLA/CARTOLA_Reaction.c \
../CARTOLA/CARTOLA_Velocity.c \
../CARTOLA/CARTOLA_upresc.c \
../HEMKER/HEMKER_Condutivity.c \
../HEMKER/HEMKER_Font.c \
../HEMKER/HEMKER_Reaction.c \
../HEMKER/HEMKER_Velocity.c \
../HEMKER/HEMKER_upresc.c \
../PUDIM/PUDIM_Condutivity.c \
../PUDIM/PUDIM_Font.c \
../PUDIM/PUDIM_Reaction.c \
../PUDIM/PUDIM_Velocity.c \
../PUDIM/PUDIM_upresc.c \
../TESTE/TESTE_Condutivity.c \
../TESTE/TESTE_Font.c \
../TESTE/TESTE_Reaction.c \
../TESTE/TESTE_Velocity.c \
../TESTE/TESTE_upresc.c \
../REACTION/REACTION_Condutivity.c \
../REACTION/REACTION_Font.c \
../REACTION/REACTION_Reaction.c \
../REACTION/REACTION_Velocity.c \
../REACTION/REACTION_upresc.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemv.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/Reordering/reordering.c \
../../../00_CommonFiles/Reordering/symrcm.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/LU_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SGS_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/SPIKE_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/Diag_precond_CSR_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/Diag_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/LU_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SGS_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/SPIKE_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/scaling.c \

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=SSTranspEquation 

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -DSSTranspEquation2D -Ddebug -c -lm -lrt -Wall -g
debug: GCFLAGS= -DSSTransEquation2D -Ddebug -lm -lrt -Wall \
../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-I/opt/openmpi/include -L/opt/openmpi/lib -lmpi -libverbs -lrt -lnsl -lutil -lm -ldl -lm -Wl,--export-dynamic \
-lrt -lnsl -lutil -lm -ldl -L/usr/local/lib -lm -lgfortran -Wall -L ../../../00_CommonFiles/Solvers_and_Preconditioners/pardiso -lpardiso500-GNU481-X86-64 \
-lblas -llapack -lquadmath -fopenmp -g 
debug: all

fast: CFLAGS= -DSSTranspEquation2D -c -lm -lrt -Wall -Ofast -march=native 
fast: GCFLAGS= -DSSTransEquation2D -lm -lrt -Wall -Ofast -march=native \
../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-I/opt/openmpi/include -L/opt/openmpi/lib -lmpi -libverbs -lrt -lnsl -lutil -lm -ldl -lm -Wl,--export-dynamic \
-lrt -lnsl -lutil -lm -ldl -L/usr/local/lib -lm -lgfortran -Wall -L ../../../00_CommonFiles/Solvers_and_Preconditioners/pardiso -lpardiso500-GNU481-X86-64 \
-lblas -llapack -lquadmath -fopenmp
fast: all

gprof: CFLAGS+=-pg -Ddebug 
gprof: LDFLAGS+=-pg -Ddebug
gprof: all
 
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)



