dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([acarpack], [1.0], [jwgcarlson@gmail.com])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE([disable])

# Check whether or not to build parallel version.
AC_ARG_ENABLE([parpack], [build Parallel ARPACK as well], [], [enable_parpack=yes])

AC_CONFIG_SUBDIRS([ARPACK])
my_subdirs="ARPACK"
if test "x$enable_parpack" = "xyes"; then
    AC_CONFIG_SUBDIRS([PARPACK])
    my_subdirs="$my_subdirs PARPACK"
fi

AC_SUBST([my_subdirs])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
