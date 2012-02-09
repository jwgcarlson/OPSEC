#!/bin/sh
#
# This script should only be run by developers, after modifying autoconf- and
# automake-related files.  It is used to generate the configure script and
# Makefile.in files needed for compilation.  Normal users should use the
# regular
#   ./configure && make
# recipe.  (Note that in principle this script should not need be run even
# after modifying build files, since the appropriate autotools calls should be
# triggered automatically by the Makefile.  It exists only for the paranoid.)

AUTORECONF=${AUTORECONF:-autoreconf}
ACLOCAL=${ACLOCAL:-aclocal}
AUTOCONF=${AUTOCONF:-autoconf}
AUTOHEADER=${AUTOHEADER:-autoheader}
AUTOMAKE=${AUTOMAKE:-automake}
LIBTOOLIZE=${LIBTOOLIZE:-libtoolize}

rm -rf autom4te.cache

if $AUTORECONF --version > /dev/null 2>&1 ; then
    echo "Using autoreconf"
    $AUTORECONF --verbose --install --symlink
else
    echo "Using manual bootstrap"
    $ACLOCAL -I m4
    $LIBTOOLIZE --install
    $AUTOCONF
    $AUTOHEADER
    $AUTOMAKE --add-missing --no-force
fi

rm -f config.cache
