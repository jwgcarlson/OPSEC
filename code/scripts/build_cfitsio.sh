#!/bin/sh
#
# Script for downloading and building CFITSIO.
#
# We assume that this is being called from OPSEC's build.sh script, so that
# all custom compilation/installation environment variables are already defined.
#
# TODO:
#  - allow for building from a separate build directory, instead of in the source tree
#  - prompt before downloading

die() {
    echo "CFITSIO build failed."
    exit 1
}

# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    die
fi

# Download tarball
if [ ! -d cfitsio ]; then
    if [ -z "$TMPDIR" ]; then
        TMPDIR=/tmp
    fi
    if [ -z "$CFITSIO_TARBALL" ]; then
        CFITSIO_TARBALL=cfitsio3290.tar.gz
    fi
    wget -P $TMPDIR -c ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/$CFITSIO_TARBALL || die
    tar -xzv -f $TMPDIR/$CFITSIO_TARBALL || die
fi


# Build
cd cfitsio || die
echo "Building all in $PWD"
./configure --prefix=$PREFIX && make && make install || die
