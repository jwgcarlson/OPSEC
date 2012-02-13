CC = $(MPICC)
CXX = $(MPICXX)

check_PROGRAMS = testMyMatrix testabn testcfg testeig testmv testsig
TESTS = $(check_PROGRAMS)

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = -L$(top_builddir)/src -lklpk

testabn_SOURCES = testabn.c CuTest.c

testcfg_SOURCES = testcfg.c CuTest.c

testeig_SOURCES = testeig.cpp CuTest.c
testeig_LDADD = $(PARPACK_LIBS) $(ARPACK_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(MPIFLIBS)

testmv_SOURCES = testmv.cpp CuTest.c
testmv_LDADD = $(ANASAZI_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)

testMyMatrix_SOURCES = testMyMatrix.cpp CuTest.c
testMyMatrix_LDADD = $(ANASAZI_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)

testsig_SOURCES = testsig.cpp
