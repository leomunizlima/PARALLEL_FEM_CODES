AC_INIT(mc75, 1.1.0, hsl@stfc.ac.uk)
AC_CONFIG_SRCDIR([src/mc75s.f])
AM_INIT_AUTOMAKE(mc75, 1.1.0)

# Check for required programs
AC_PROG_F77
AC_PROG_FC
AC_PROG_CC
AC_PROG_RANLIB

# Establish some required flags
AC_FC_LIBRARY_LDFLAGS # Linking C to Fortran

# Establish correct module search path flag
AX_F90_MODULE_FLAG

# Check for external libraries (e.g. BLAS, LAPACK, MeTiS)
ACX_BLAS(,[AC_MSG_WARN([No BLAS library found. make check will not work.])])

# Output data
AC_CONFIG_FILES(Makefile src/Makefile example/Makefile)
AC_OUTPUT
