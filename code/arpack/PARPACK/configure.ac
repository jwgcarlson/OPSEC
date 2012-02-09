dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT([PARPACK], [2.1], [jwgcarlson@gmail.com])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE([disable])
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_LIBTOOL
AC_PROG_F77

# Find MPI-compatible F77 compiler
AC_LANG([Fortran 77])
AX_MPI

F77="$MPIF77"
LIBS="$MPILIBS $LIBS"

AC_CONFIG_FILES([Makefile SRC/MPI/Makefile])
AC_OUTPUT
