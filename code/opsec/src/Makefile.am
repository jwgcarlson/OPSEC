CC = $(MPICC)
CXX = $(MPICXX)
F77 = $(MPIF77)

# Static convenience library.
noinst_LTLIBRARIES = libklpk.la

libklpk_la_SOURCES = \
    MyMatrix.hpp MyMultiVec.h MySignalMatrixOperator.h \
    vec3.h Matrix.h \
    abn.c cfg.c chealpix.c eig.c opsec.c rng.c sig.cpp Cell.cpp QSpline.cpp \
    abn.h cfg.h chealpix.h eig.h opsec.h rng.h sig.h   Cell.h   QSpline.h \
    Spline.cpp Survey.cpp BoxSurvey.cpp SphericalSurvey.cpp \
    Spline.h Survey.h   BoxSurvey.h   SphericalSurvey.h \
    Model.cpp KaiserModel.cpp RealModel.cpp RealSplineModel.cpp SzapudiModel.cpp TestModel.cpp XiModel.cpp \
    Model.h   KaiserModel.h   RealModel.h   RealSplineModel.h   SzapudiModel.h   TestModel.h   XiModel.h
libklpk_la_CXXFLAGS = $(OPENMP_CFLAGS)

bin_PROGRAMS = basis comma dot estimate klt klt-anasazi

basis_SOURCES = basis.cpp
basis_CXXFLAGS = $(OPENMP_CFLAGS)
basis_LDADD = libklpk.la $(FITSIO_LIBS)

comma_SOURCES = comma.cpp
comma_CXXFLAGS = $(OPENMP_CFLAGS)
comma_LDADD = libklpk.la $(BLAS_LIBS) $(FLIBS)

dot_SOURCES = dot.cpp
dot_CXXFLAGS = $(OPENMP_CFLAGS)
dot_LDADD = libklpk.la $(FITSIO_LIBS) $(BLAS_LIBS) $(FLIBS)

estimate_SOURCES = estimate.cpp
estimate_CXXFLAGS = $(OPENMP_CFLAGS)
estimate_LDADD = libklpk.la $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

klt_SOURCES = klt.cpp
klt_CXXFLAGS = $(OPENMP_CFLAGS)
klt_LDADD = libklpk.la $(PARPACK_LIBS) $(ARPACK_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(MPIFLIBS)

klt_anasazi_SOURCES = klt-anasazi.cpp
klt_anasazi_CXXFLAGS = $(OPENMP_CFLAGS)
klt_anasazi_LDADD = libklpk.la $(ANASAZI_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)


printflags:
	@echo FITSIO_LIBS = $(FITSIO_LIBS)
	@echo ANASAZI_LIBS = $(ANASAZI_LIBS)
	@echo ARPACK_LIBS = $(ARPACK_LIBS)
	@echo PARPACK_LIBS = $(PARPACK_LIBS)
	@echo BLAS_LIBS = $(BLAS_LIBS)
	@echo LAPACK_LIBS = $(LAPACK_LIBS)
	@echo FLIBS = $(FLIBS)
	@echo MPIFLIBS = $(MPIFLIBS)
