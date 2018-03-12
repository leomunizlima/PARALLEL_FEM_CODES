CC=mpicc
#CC=tau_cc.sh
IDIR=
CFLAGS= -DNavierStokesEquations2D -Ddebug -c -lm -lrt -Wall -Ofast -march=native
GCFLAGS= -DNavierStokesEquations2D -Ddebug -lm -lrt -Wall -Ofast -march=native ../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
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
Fill_LM.c \
setProblem.c \
setSolver.c \
setStabilizationForm.c \
setPreconditioner.c \
setMatrixVectorProductType.c \
setzeros.c \
ebe_assembly.c \
csr_Initialization.c \
csr_assembly.c \
COMPARE_array_in_Preprocess.c \
MPI_VectorUpdate.c \
eval_U.c \
Build_K_F_SUPG_PSPG.c \
Build_K_F_VMS_DCDD.c \
../CAVITY/CAVITY_v1presc.c \
../CAVITY/CAVITY_v2presc.c \
../CAVITY/CAVITY_ppresc.c \
../CHANNEL/CHANNEL_v1presc.c \
../CHANNEL/CHANNEL_v2presc.c \
../CHANNEL/CHANNEL_ppresc.c \
../EXATA/EXATA_v1presc.c \
../EXATA/EXATA_v2presc.c \
../EXATA/EXATA_ppresc.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/Diag_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/Diag_precond_CSR_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/BlockDiagDOF3_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/SPIKE_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/BlockDiagDOF3_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/SPIKE_precond.c \
../../../00_CommonFiles/BLAS_Operations/dtrsvUP.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/BLAS_Operations/izero.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemvNDOF3.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/Reordering/reordering.c \
../../../00_CommonFiles/Reordering/symrcm.c \

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=NavierStokesEquations2D

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -DNavierStokesEquations2D -g -Ddebug -c -lm -lrt -Wall
debug: GCFLAGS= -DNavierStokesEquations2D -g -Ddebug -lm -lrt -Wall ../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-I/opt/openmpi/include -L/opt/openmpi/lib -lmpi -libverbs -Wl,--export-dynamic \
-lnsl -lutil -ldl -L/usr/local/lib -lgfortran -L ../../../00_CommonFiles/Solvers_and_Preconditioners/pardiso -lpardiso500-GNU481-X86-64 \
-lblas -llapack -lquadmath -fopenmp
debug: all

fast: CFLAGS= -DNavierStokesEquations2D -c -lm -lrt -Wall -Ofast -march=native
fast: GCFLAGS= -DNavierStokesEquations2D -lm -lrt -Wall -Ofast -march=native ../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-I/opt/openmpi/include -L/opt/openmpi/lib -lmpi -libverbs -Wl,--export-dynamic \
-lnsl -lutil -ldl -L/usr/local/lib -lgfortran -L ../../../00_CommonFiles/Solvers_and_Preconditioners/pardiso -lpardiso500-GNU481-X86-64 \
-lblas -llapack -lquadmath -fopenmp
fast: all

gprof: CFLAGS+=-pg
gprof: LDFLAGS+=-pg
gprof: all


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
