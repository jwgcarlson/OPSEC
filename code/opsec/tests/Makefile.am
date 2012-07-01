CC = $(MPICC)
CXX = $(MPICXX)

check_PROGRAMS = testBoxSurvey testMyMatrix testSplitFile testSurvey \
                 testabn testcfg testeig testmv testsig testslp
TESTS = $(check_PROGRAMS)

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = -L$(top_builddir)/src -lklpk

testabn_SOURCES = testabn.c CuTest.c

testcfg_SOURCES = testcfg.c CuTest.c

testeig_SOURCES = testeig.cpp CuTest.c
testeig_LDADD = $(PARPACK_LIBS) $(ARPACK_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(MPIFLIBS)

testmv_SOURCES = testmv.cpp CuTest.c
testmv_LDADD = $(ANASAZI_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)

testsig_SOURCES = testsig.cpp CuTest.c

testBoxSurvey_SOURCES = testBoxSurvey.cpp CuTest.c

testMyMatrix_SOURCES = testMyMatrix.cpp CuTest.c
testMyMatrix_LDADD = $(ANASAZI_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)

testSplitFile_SOURCES = testSplitFile.cpp CuTest.c

testSurvey_SOURCES = testSurvey.cpp CuTest.c

testslp_SOURCES = testslp.cpp CuTest.c
testslp_LDADD = $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)
