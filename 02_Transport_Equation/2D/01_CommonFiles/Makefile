CC=mpicc
#CC=tau_cc.sh
IDIR= 
CFLAGS= -DTranspEquation2D -Ddebug -c -lm -lrt -Wall -Ofast -march=native 
GCFLAGS= -DTransEquation2D -Ddebug -lm -lrt -Wall -Ofast -march=native \
../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-I/opt/openmpi/include -L/opt/openmpi/lib -lmpi -libverbs -Wl,--export-dynamic \
-lnsl -lutil -ldl -L/usr/local/lib -lgfortran -L ../../../00_CommonFiles/Solvers_and_Preconditioners/pardiso -lpardiso500-GNU481-X86-64 \
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
setPreconditioner.c \
setStabilizationForm.c \
setzeros.c \
csr_Initialization.c \
ebe_assembly.c \
csr_assembly.c \
F_assembly.c \
eval_U_dU.c \
Build_M_K_R_SUPG.c \
Build_M_K_R_DD.c \
Build_M_F_DD.c \
calculate_DaB.c \
CAU_ShockCapture.c \
CAU_DD_ShockCapture.c \
YZBeta_ShockCapture.c \
Fill_LM.c \
Permutation_of_LM.c \
h_shock_2sqrtArea.c \
h_shock_Option1.c \
h_shock_Option2.c \
MPI_VectorUpdate.c \
COMPARE_array_in_Preprocess.c \
../CARTOLA/CARTOLA_Condutivity.c \
../CARTOLA/CARTOLA_Font.c \
../CARTOLA/CARTOLA_Reaction.c \
../CARTOLA/CARTOLA_Velocity.c \
../CARTOLA/CARTOLA_upresc.c \
../CARTOLA/CARTOLA_InitialSolution.c \
../CONE/CONE_Condutivity.c \
../CONE/CONE_Font.c \
../CONE/CONE_Reaction.c \
../CONE/CONE_Velocity.c \
../CONE/CONE_upresc.c \
../CONE/CONE_InitialSolution.c \
../CONE2/CONE2_Condutivity.c \
../CONE2/CONE2_Font.c \
../CONE2/CONE2_Reaction.c \
../CONE2/CONE2_Velocity.c \
../CONE2/CONE2_upresc.c \
../CONE2/CONE2_InitialSolution.c \
../PUDIM/PUDIM_Condutivity.c \
../PUDIM/PUDIM_Font.c \
../PUDIM/PUDIM_Reaction.c \
../PUDIM/PUDIM_Velocity.c \
../PUDIM/PUDIM_upresc.c \
../PUDIM/PUDIM_InitialSolution.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/C_Operations/dmemcpy.c \
../../../00_CommonFiles/C_Operations/memsetzero.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemv.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/Reordering/reordering.c \
../../../00_CommonFiles/Reordering/symrcm.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/SPIKE_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/Diag_precond_CSR_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/Diag_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/SPIKE_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
../../../00_CommonFiles/Time_Integration/Criteria.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old.c \
../../../00_CommonFiles/Time_Integration/Predictor_New.c \

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=TranspEquation 

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -DTranspEquation2D -Ddebug -c -lm -lrt -Wall -g
debug: GCFLAGS= -DTransEquation2D -Ddebug -lm -lrt -Wall \
../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-I/opt/openmpi/include -L/opt/openmpi/lib -lmpi -libverbs -Wl,--export-dynamic \
-lnsl -lutil -ldl -L/usr/local/lib -lgfortran -L ../../../00_CommonFiles/Solvers_and_Preconditioners/pardiso -lpardiso500-GNU481-X86-64 \
-lblas -llapack -lquadmath -fopenmp -g 
debug: all

fast: CFLAGS= -DTranspEquation2D -c -lm -lrt -Wall -Ofast -march=native 
fast: GCFLAGS= -DTransEquation2D -lm -lrt -Wall -Ofast -march=native \
../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-I/opt/openmpi/include -L/opt/openmpi/lib -lmpi -libverbs -Wl,--export-dynamic \
-lnsl -lutil -ldl -L/usr/local/lib -lgfortran -L ../../../00_CommonFiles/Solvers_and_Preconditioners/pardiso -lpardiso500-GNU481-X86-64 \
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



