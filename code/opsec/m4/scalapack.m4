# Detect ScaLAPACK library.
#
# Original version taken from the Octopus project
# (http://www.tddft.org/programs/octopus), version 4.0.1.
#
# Modified by Jordan Carlson <jwgcarlson@gmail.com>, 18-2-2012.

AC_DEFUN([ACX_SCALAPACK], [
AC_REQUIRE([ACX_BLACS])
acx_scalapack_ok=no

# We cannot use ScaLAPACK if BLACS is not found.
if test "x$acx_blacs_ok" != xyes; then
    acx_scalapack_ok=noblacs
fi

# Get Fortran linker name of ScaLAPACK function to check for
AC_F77_FUNC(pcheev)

# Check if ScaLAPACK library was given explicitly
if test $acx_scalapack_ok = no; then
    AC_ARG_WITH(scalapack, [AS_HELP_STRING([--with-scalapack=<lib>], [use SCALAPACK library <lib>])])
    case $with_scalapack in
        yes | "") ;;
        no) acx_scalapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LIBS_SCALAPACK="$with_scalapack" ;;
        *) LIBS_SCALAPACK="-l$with_scalapack" ;;
    esac
fi

# Backup LIBS 
acx_scalapack_save_LIBS="$LIBS"
LIBS="$BLACS_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

# First, check SCALAPACK_LIBS environment variable
if test $acx_scalapack_ok = no; then
    if test "x$SCALAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$SCALAPACK_LIBS $LIBS"
        AC_MSG_CHECKING([for $pcheev in $SCALAPACK_LIBS])
        AC_TRY_LINK_FUNC($pcheev, [acx_scalapack_ok=yes], [SCALAPACK_LIBS=""])
        AC_MSG_RESULT([$acx_scalapack_ok ($SCALAPACK_LIBS)])
        LIBS="$save_LIBS"
    fi
fi

# ScaLAPACK linked to by default?  (for instance if given in BLAS_LIBS)
if test $acx_scalapack_ok = no; then
    AC_CHECK_FUNC($pcheev, [acx_scalapack_ok=yes])
fi

# Generic ScaLAPACK library?
for scalapack in scalapack scalapack-openmpi; do
    if test $acx_scalapack_ok = no; then
        AC_CHECK_LIB($scalapack, $pcheev,
            [acx_scalapack_ok=yes; SCALAPACK_LIBS="$SCALAPACK_LIBS -l$scalapack"], [], [$FLIBS])
    fi
done

AC_SUBST(SCALAPACK_LIBS)
LIBS="$acx_scalapack_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_scalapack_ok" = xyes; then
    ifelse([$1],,AC_DEFINE(HAVE_SCALAPACK,1,[Define if you have ScaLAPACK library.]),[$1])
else
    acx_scalapack_ok=no
    $2
fi
])dnl ACX_SCALAPACK
