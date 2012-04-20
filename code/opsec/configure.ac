#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([opsec], [0.1], [jwgcarlson@gmail.com])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE([disable])
LT_INIT([disable-shared])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([opsec_config.h])
AC_CONFIG_SRCDIR([src/cfg.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_WRAPPERS
AC_PROG_LIBTOOL

# Check for Python.
AS_IF([test -z "$PYTHON" -a -n "$OPSEC_PYTHON"],
      [PYTHON=`which "$OPSEC_PYTHON"`]
)
AM_PATH_PYTHON([2.5], [], [AC_MSG_WARN([Python support disabled.])])
AM_CONDITIONAL([OPSEC_AM_HAVE_PYTHON], [test "$PYTHON" != ":"])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr strdup strndup])

# Debug build?
AC_ARG_ENABLE([debug],
              [--enable-debug   Turn on debugging],
              [case "$enableval" in
                   yes) debug=true ;;
                   no)  debug=false ;;
                   *) AC_MSG_ERROR([bad value $enableval for --enable-debug]) ;;
               esac],
              [debug=false])
AS_IF([test x$debug = xtrue],
      [AC_DEFINE([OPSEC_DEBUG], 1, [Define to enable debugging])]
)
AM_CONDITIONAL([OPSEC_AM_DEBUG], [test x$debug = xtrue])

# Decide on precision: single or double.
AC_ARG_WITH(
    [precision],
    [AS_HELP_STRING([--with-precision=PREC], [choose single or double (default) precision])],
    [PRECISION="$withval"],
    [PRECISION="double"]
)
if test "$PRECISION" = "single"; then
    AC_DEFINE(OPSEC_SINGLE, 1, [Use single precision floating point operations])
elif test "$PRECISION" = "double"; then
    AC_DEFINE(OPSEC_DOUBLE, 1, [Use double precision floating point operations])
else
    AC_MSG_ERROR([precision must be "single" or "double", not "$PRECISION"])
fi


dnl *******************************************************************
dnl Find MPI-compatible compilers
dnl *******************************************************************

# OPSEC currently doesn't compile without MPI support, and it wouldn't be very
# useful even if it did.  So let's stop pretending to be so flexible and just
# assume that MPI is present.
AC_DEFINE(OPSEC_USE_MPI, 1, [Define if OPSEC should use MPI.])

# In an ideal world, we'd just use whatever MPI-compatible compilers are
# available on our system.
AC_LANG_PUSH([C])
    AX_MPI([have_mpicc=yes],
           [AC_MSG_ERROR([cannot compile MPI program in C, try setting MPICC])])
AC_LANG_POP([C])
AC_LANG_PUSH([C++])
    AX_MPI([have_mpicxx=yes],
           [AC_MSG_ERROR([cannot compile MPI program in C++, try setting MPICXX])])
AC_LANG_POP([C++])
AC_LANG_PUSH([Fortran 77])
    AX_MPI([have_mpif77=yes],
           [AC_MSG_ERROR([cannot compile MPI program in F77, try setting MPIF77])])
AC_LANG_POP([Fortran 77])

# Unfortunately, we need to link C++ programs with Fortran libraries.  With
# OpenMPI, this requires you to use mpicxx as your linker, with the additional
# linker flag -lmpi_f77.  I'm not sure what magic needs to be done with other
# compilers, but hopefully it's not too much harder.
AC_ARG_VAR(MPIFLIBS, [additional linker flags to link MPI-enabled Fortran libraries in C/C++ programs])
if test "x$MPIFLIBS" = x; then
    case "$MPICXX" in
    *mpicxx*|*mpic++*|*mpiCC*)
        AC_LANG_PUSH([C++])
            save_CXX="$CXX"
            save_LIBS="$LIBS"
            CXX="$MPICXX"
            LIBS="$FLIBS $MPILIBS $LIBS"
            AC_CHECK_LIB([mpi_f77], [MPI_INIT], [MPIFLIBS="-lmpi_f77"])
            CXX="$save_CXX"
            LIBS="$save_LIBS"
        AC_LANG_POP([C++])
        ;;
    esac
fi
AC_SUBST(MPIFLIBS)


dnl *******************************************************************
dnl Check for OpenMP.
dnl *******************************************************************

AC_OPENMP


dnl *******************************************************************
dnl Find BLAS and LAPACK
dnl *******************************************************************

AX_BLAS([], AC_MSG_ERROR([Could not find BLAS.]))
AX_LAPACK([], AC_MSG_ERROR([Could not find LAPACK.]))


dnl *******************************************************************
dnl Find BLACS and ScaLAPACK
dnl *******************************************************************

AC_LANG_PUSH([C])
    save_CC="$CC"
    CC="$MPICC"
    ACX_BLACS([], AC_MSG_ERROR([Could not find BLACS.]))
    ACX_SCALAPACK([], AC_MSG_ERROR([Could not find ScaLAPACK.]))
    CC="$save_CC"
AC_LANG_POP([C])


dnl *******************************************************************
dnl Find CFITSIO headers and library
dnl *******************************************************************

AC_ARG_WITH(
    [cfitsio-include],
    [AS_HELP_STRING([--with-cfitsio-include=PATH], [path to CFITSIO headers])],
    [FITSIO_CPPFLAGS="-I$withval"],
    [FITSIO_CPPFLAGS=""]
)
AC_ARG_WITH(
    [cfitsio-libs],
    [AS_HELP_STRING([--with-cfitsio-libs=PATH], [linker options for CFITSIO])],
    [FITSIO_LIBS="$withval"],
    [FITSIO_LIBS="-lcfitsio"]
)
save_CPPFLAGS="$CPPFLAGS"
save_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $FITSIO_CPPFLAGS"
LIBS="$LIBS $FITSIO_LIBS -lm"
AC_CHECK_HEADER([fitsio.h],
                [],
                [AC_MSG_FAILURE([Could not find CFITISO headers. Try setting --with-fitsio-include.])]
)
AC_TRY_LINK([#include <fitsio.h>],
            [fits_init_cfitsio();],
            [],
            [AC_MSG_FAILURE([Could not link against CFITSIO library. Try setting --with-fitsio-libs.])]
)
CPPFLAGS="$save_CPPFLAGS"
LIBS="$save_LIBS"
AC_SUBST(FITSIO_CPPFLAGS)
AC_SUBST(FITSIO_LIBS)


dnl *******************************************************************
dnl Find ARPACK/PARPACK libraries
dnl *******************************************************************

AC_ARG_WITH(
    [arpack],
    [AS_HELP_STRING([--with-arpack=PATH], [full path to libarpack.a])],
    [ARPACK_LIBS="$with_arpack"],
    [ARPACK_LIBS="-larpack"]
)
AC_ARG_WITH(
    [parpack],
    [AS_HELP_STRING([--with-parpack=PATH], [full path to libparpack.a])],
    [PARPACK_LIBS="$with_parpack"],
    [PARPACK_LIBS="-lparpack"]
)

AC_LANG_PUSH([Fortran 77])
save_LIBS="$LIBS"
LIBS="$ARPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $FLIBS $LIBS"
AC_LINK_IFELSE(AC_LANG_CALL([], [dseupd]),
               [AC_DEFINE(HAVE_ARPACK, 1, [Define if ARPACK is available.])],
               [AC_MSG_FAILURE([Could not link against ARPACK library. Try setting --with-arpack.])]
)
save_F77="$F77"
F77="$MPIF77"
LIBS="$PARPACK_LIBS $LIBS"
AC_LINK_IFELSE(AC_LANG_CALL([], [pdseupd]),
               [AC_DEFINE(HAVE_PARPACK, 1, [Define if PARPACK is available.])],
               [AC_MSG_FAILURE([Could not link against PARPACK library. Try setting --with-parpack.])]
)
LIBS="$save_LIBS"
F77="$save_F77"
AC_LANG_POP([Fortran 77])

AC_SUBST(ARPACK_LIBS)
AC_SUBST(PARPACK_LIBS)


dnl *******************************************************************
dnl Find Anasazi.
dnl *******************************************************************

AC_ARG_WITH(
    [anasazi-include],
    [AS_HELP_STRING([--with-anasazi-include=PATH], [path to Anasazi headers])],
    [ANASAZI_CPPFLAGS="-I$withval"],
    [ANASAZI_CPPFLAGS=""]
)
AC_ARG_WITH(
    [anasazi-libs],
    [AS_HELP_STRING([--with-anasazi-libs=PATH], [linker options for Anasazi library])],
    [ANASAZI_LIBS="$withval"],
    [ANASAZI_LIBS="-lanasazi -lteuchos"]
)
save_CPPFLAGS="$CPPFLAGS"
save_CXX="$CXX"
save_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $ANASAZI_CPPFLAGS"
CXX="$MPICXX"
LIBS="$LIBS $ANASAZI_LIBS -lm"
AC_LANG_PUSH([C++])
AC_CHECK_HEADER([AnasaziBasicSort.hpp],
                [],
                [AC_MSG_FAILURE([Could not find Anasazi headers. Try setting --with-anasazi-include.])],
                AC_INCLUDES_DEFAULT
)
AC_MSG_CHECKING([Anasazi usability])
AC_LINK_IFELSE(
    AC_LANG_PROGRAM(
        [[#include <AnasaziBasicSort.hpp>
          #include <Teuchos_ParameterList.hpp>]],
        [[Teuchos::ParameterList pl;
          Anasazi::BasicSort<double> sort(pl);]]
    ),
    [AC_MSG_RESULT(yes)
     AC_DEFINE(HAVE_ANASAZI, 1, [Define if Anasazi is available.])],
    [AC_MSG_RESULT(no)
     AC_MSG_FAILURE([Could not link against Anasazi library. Try setting --with-anasazi-libs.])]
)
AC_LANG_POP([C++])
CPPFLAGS="$save_CPPFLAGS"
CXX="$save_CXX"
LIBS="$save_LIBS"
AC_SUBST(ANASAZI_CPPFLAGS)
AC_SUBST(ANASAZI_LIBS)


dnl *******************************************************************
dnl Check for CImg.
dnl *******************************************************************

AC_LANG_PUSH([C++])
AC_CHECK_HEADER(CImg.h, cimg_ok=yes, cimg_ok=no)
AC_LANG_POP([C++])
AM_CONDITIONAL([OPSEC_AM_HAVE_CIMG], [test x$cimg_ok = xyes])


dnl *******************************************************************
dnl Check for libpng.
dnl *******************************************************************

AC_CHECK_LIB(png, png_read_info,
    AC_CHECK_HEADER(png.h, png_ok=yes, png_ok=no),
    png_ok=no,
    -lz
)
AM_CONDITIONAL([OPSEC_AM_HAVE_LIBPNG], [test x$png_ok = xyes])
AC_SUBST(LIBPNG, ["-lpng -lz"])


dnl *******************************************************************
dnl Find FFTW2 headers and libraries, for gmock utility program
dnl *******************************************************************

FFTW2_LIBS="-lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -lm"
save_CC="$CC"
save_LIBS="$LIBS"
CC="$MPICC"
LIBS="$LIBS $FFTW2_LIBS"
AC_LINK_IFELSE(
    AC_LANG_PROGRAM([[#include <rfftw.h>
                      #include <rfftw_mpi.h>]],
                    [[rfftwnd_mpi_destroy_plan(NULL);]]),
    fftw2_ok=yes,
    fftw2_ok=no
)
CC="save_CC"
LIBS="$save_LIBS"

AM_CONDITIONAL([OPSEC_AM_HAVE_FFTW2], [test x$fftw2_ok = xyes])
AC_SUBST(FFTW2_LIBS)


dnl *******************************************************************
dnl Write Makefiles
dnl *******************************************************************

AC_SUBST(CPPFLAGS, ["$CPPFLAGS $ANASAZI_CPPFLAGS $FITSIO_CPPFLAGS"])
AC_SUBST(UTILITIES)

AC_CONFIG_FILES([Makefile
                 python/Makefile
                 python/opsec-submit
                 src/Makefile
                 tests/Makefile
                 utils/Makefile])
AC_OUTPUT


dnl *******************************************************************
dnl Summarize configuration
dnl *******************************************************************
cat << EOF
------------------------------------------------------------------------
OPSEC configuration summary:

CC = $CC
CXX = $CXX
F77 = $F77
MPICC = $MPICC
MPICXX = $MPICXX
MPIF77 = $MPIF77
CPPFLAGS = $CPPFLAGS
CFLAGS = $CFLAGS
CXXFLAGS = $CXXFLAGS
F77FLAGS = $F77FLAGS
LDFLAGS = $LDFLAGS
LIBS = $LIBS
FLIBS = $FLIBS
MPILIBS = $MPILIBS
MPIFLIBS = $MPIFLIBS
UTILITIES = $UTILITIES
------------------------------------------------------------------------
EOF
