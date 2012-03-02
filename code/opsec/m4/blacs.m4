# Detect BLACS library.
#
# Original version taken from the Octopus project
# (http://www.tddft.org/programs/octopus), version 4.0.1.
#
# Modified by Jordan Carlson <jwgcarlson@gmail.com>, 18-2-2012.

AC_DEFUN([ACX_BLACS], [
AC_REQUIRE([AX_MPI])
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_blacs_ok=no

# Check if BLACS library was given explicitly
AC_ARG_WITH(blacs, [AS_HELP_STRING([--with-blacs=<lib>], [use BLACS library <lib>])])
case $with_blacs in
    yes | "") ;;
    no) acx_blacs_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) BLACS_LIBS="$with_blacs" ;;
    *) BLACS_LIBS="-l$with_blacs" ;;
esac

# Get Fortran linker name of BLACS function to check for
AC_F77_FUNC(blacs_pinfo)

# Backup LIBS 
acx_blacs_save_LIBS="$LIBS"
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_MSG_NOTICE([blacs.m4: BLACS_LIBS = $BLACS_LIBS])
# First, check BLACS_LIBS environment variable
if test $acx_blacs_ok = no; then
    if test "x$BLACS_LIBS" != x; then
        AC_MSG_NOTICE([blacs.m4: yep, BLACS_LIBS = $BLACS_LIBS])
        save_LIBS="$LIBS"; LIBS="$BLACS_LIBS $LIBS"
        AC_MSG_CHECKING([for $blacs_pinfo in $BLACS_LIBS])
        AC_TRY_LINK_FUNC($blacs_pinfo, [acx_blacs_ok=yes], [BLACS_LIBS=""])
        AC_MSG_RESULT([$acx_blacs_ok ($BLACS_LIBS)])
        LIBS="$save_LIBS"
    fi
fi

# BLACS linked to by default?  (for instance if given in BLAS_LIBS)
if test $acx_blacs_ok = no; then
    AC_CHECK_FUNC($blacs_pinfo, [acx_blacs_ok=yes])
fi

# Generic BLACS library?
for blacs in blacs blacs-openmpi; do
    if test x"$blacs" = xblacs-openmpi; then       
        blacsinit="blacsF77init-openmpi"
    else
        blacsinit="blacsF77init"
    fi
    if test $acx_blacs_ok = no; then
        AC_CHECK_LIB($blacs -l$blacsinit -l$blacs, $blacs_pinfo,
            [acx_blacs_ok=yes; BLACS_LIBS="$BLACS_LIBS -l$blacs -l$blacsinit -l$blacs"], [], [$FLIBS])
    fi
done

AC_SUBST(BLACS_LIBS)
LIBS="$acx_blacs_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blacs_ok" = xyes; then
    ifelse([$1],,AC_DEFINE(HAVE_BLACS,1,[Define if you have BLACS library.]),[$1])
else
    acx_blacs_ok=no
    $2
fi
])dnl ACX_BLACS
