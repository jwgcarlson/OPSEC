#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([gmock], [0.1], [jwgcarlson@gmail.com])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE([disable])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/cfg.h])

# See if MPI was requested.
AC_ARG_ENABLE([mpi],
              [AS_HELP_STRING([--enable-mpi], [compile with MPI support (default yes)])],
              [USE_MPI="$enableval"],
              [USE_MPI=yes])

# Allow user to manually override compiler settings.
if test -z "$CC" && test -n "$MPICC"; then
    CC="$MPICC"
fi
if test -z "$CXX" && test -n "$MPICXX"; then
    CXX="$MPICXX"
fi

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL

# Check for MPI availability.
if test x"$USE_MPI" = xyes; then
    AC_LANG_PUSH([C])
        AX_MPI([], AC_MSG_ERROR([could not compile MPI test program in C]))
    AC_LANG_POP([C])
    AC_LANG_PUSH([C++])
        AX_MPI([], AC_MSG_ERROR([could not compile MPI test program in C++]))
    AC_LANG_POP([C++])

    CC="$MPICC"
    CXX="$MPICXX"
    LIBS="$MPILIBS $LIBS"
fi

# Checks for generic libraries.

# Checks for header files.

# Checks for FFTW 2.
AC_CHECK_HEADER([rfftw.h], [], [AC_MSG_ERROR([could not find <rfftw.h>, try setting CPPFLAGS])])
AC_CHECK_LIB([fftw], [fftw_malloc], [], [AC_MSG_ERROR([could not link -lfftw, try setting LDFLAGS])])
AC_CHECK_LIB([rfftw], [rfftw], [], [AC_MSG_ERROR([could not link -lrfftw, try setting LDFLAGS])])
if test x"$USE_MPI" = xyes; then
    AC_CHECK_HEADER([rfftw_mpi.h], [], [AC_MSG_ERROR([could not find <rfftw_mpi.h>, try setting CPPFLAGS])])
    AC_CHECK_LIB([fftw_mpi], [fftw_mpi], [], [AC_MSG_ERROR([could not link -lfftw_mpi, try setting LDFLAGS])])
    AC_CHECK_LIB([rfftw_mpi], [rfftwnd_mpi], [], [AC_MSG_ERROR([could not link -lrfftw_mpi, try setting LDFLAGS])])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor memset pow sqrt strchr strdup strndup])
AC_FUNC_MALLOC
AC_FUNC_STRTOD

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
