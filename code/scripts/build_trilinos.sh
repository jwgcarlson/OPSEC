#!/bin/sh
#
# Script for downloading and building Trilinos.
#
# We assume that this is being called from OPSEC's build.sh script, so that
# all custom compilation/installation environment variables are already defined.
#
# TODO:
#  - allow for building from a separate build directory, instead of in the source tree
#  - prompt before downloading


# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    die
fi

die() {
    echo "Trilinos build failed."
    exit 1
}

# Ensure CMake is available on the system
if [ -z "$CMAKE" ]; then
    CMAKE=cmake
fi
if [ -n "`which $CMAKE`" ] && $CMAKE --version ; then
    :
else
    echo "CMake not found.  Please install it (www.cmake.org) and try again."
    echo "If installed to a non-standard path, set the environment variable CMAKE."
    die
fi

# Choose Trilinos version
if [ -z "$TRILINOS_VERSION" ]; then
    TRILINOS_VERSION=10.8.5
fi
TRILINOS_SRCDIR="trilinos-$TRILINOS_VERSION-Source"

# Download tarball
if [ ! -d "$TRILINOS_SRCDIR" ]; then
    if [ -z "$TMPDIR" ]; then
        TMPDIR=/tmp
    fi
    TRILINOS_TARBALL=trilinos-$TRILINOS_VERSION-Source.tar.gz
    wget -P $TMPDIR -c http://trilinos.sandia.gov/download/files/$TRILINOS_TARBALL || die
    tar -xzv -f $TMPDIR/$TRILINOS_TARBALL || die
fi


# Build
cd $TRILINOS_SRCDIR || die
echo "Building all in $PWD"
mkdir -p build
cd build
if [-e CMakeCache.txt ]; then
    echo "Removing CMakeCache.txt"
    rm -f CMakeCache.txt
fi
$CMAKE \
    -D CMAKE_INSTALL_PREFIX:PATH="$PREFIX" \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D TPL_ENABLE_Boost:BOOL=ON \
    -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=FALSE \
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=FALSE \
    -D Trilinos_ENABLE_Anasazi:BOOL=ON \
    -D Trilinos_ENABLE_Teuchos:BOOL=ON \
    ../ || die
make && make install || die
