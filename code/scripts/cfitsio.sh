# Shell functions for downloading/building/cleaning CFITSIO source code.
#
# We assume these methods are being called from build.sh or clean.sh, so that
# setup.sh has been sourced and all relevant configuration variables defined.

# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    exit 1
fi

# Choose cfitsio version (default to 3290)
if [ -z "$CFITSIO_VERSION" ]; then
    CFITSIO_VERSION=3290
fi
if [ -z "$CFITSIO_DOWNLOAD_SITE" ]; then
    CFITSIO_DOWNLOAD_SITE="ftp://heasarc.gsfc.nasa.gov/software/fitsio/c"
fi
if [ -z "$CFITSIO_SRCDIR" ]; then
    CFITSIO_SRCDIR="cfitsio"
fi

download_cfitsio() {
    if [ -z "$TMPDIR" ]; then
        TMPDIR=/tmp
    fi
    CFITSIO_TARBALL=cfitsio$CFITSIO_VERSION.tar.gz
    wget -P $TMPDIR -c $CFITSIO_DOWNLOAD_SITE/$CFITSIO_TARBALL || exit 1
    tar -xzv -f $TMPDIR/$CFITSIO_TARBALL || exit 1
}

build_cfitsio() {
    if [ ! -d "$CFITSIO_SRCDIR" ]; then
        download_cfitsio
    fi
    cd $CFITSIO_SRCDIR || exit 1
    ./configure --prefix="$PREFIX" || exit 1
    make || exit 1
    make install || exit 1
}

clean_cfitsio() {
    target=clean
    if [ -n "$1" ]; then
        target="$1"
    fi
    cd $CFITSIO_SRCDIR || exit 1
    if [ -e Makefile ]; then
        make $target || exit 1
    fi
}
