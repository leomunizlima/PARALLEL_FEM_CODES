AM_FCFLAGS = -I $(top_builddir)/src

TESTS_ENVIRONMENT = $(srcdir)/runtest.sh
EXTRA_DIST = runtest.sh
TESTS = \
	mc60ss \
	mc60ds

check_PROGRAMS = \
	mc60ss \
	mc60ds

mc60ss_SOURCES = mc60ss.f
mc60ss_LDADD = $(top_builddir)/src/libmc60.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += mc60ss.output mc60ss.data
mc60ds_SOURCES = mc60ds.f
mc60ds_LDADD = $(top_builddir)/src/libmc60.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += mc60ds.output mc60ds.data
