AM_FCFLAGS = -I $(top_builddir)/src

TESTS_ENVIRONMENT = $(srcdir)/runtest.sh
EXTRA_DIST = runtest.sh
TESTS = \
	mc75ss \
	mc75ds

check_PROGRAMS = \
	mc75ss \
	mc75ds

mc75ss_SOURCES = mc75ss.f
mc75ss_LDADD = $(top_builddir)/src/libmc75.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += mc75ss.output mc75ss.data
mc75ds_SOURCES = mc75ds.f
mc75ds_LDADD = $(top_builddir)/src/libmc75.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += mc75ds.output mc75ds.data
