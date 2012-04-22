# Shell functions for downloading/building/cleaning Trilinos source code.
#
# We assume these methods are being called from build.sh or clean.sh, so that
# setup.sh has been sourced and all relevant configuration variables defined.

# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    exit 1
fi

# Ensure CMake is available on the system
if [ -z "$CMAKE" ]; then
    CMAKE=cmake
fi
if [ -n "`which $CMAKE`" ] && $CMAKE --version ; then
    :
else
    echo "CMake not found.  Please install it (www.cmake.org) and try again."
    echo "If installed to a non-standard path, set the environment variable CMAKE."
    exit 1
fi

# Choose Trilinos version (default to 10.10.1)
if [ -z "$TRILINOS_VERSION" ]; then
    TRILINOS_VERSION=10.10.1
fi

TRILINOS_SRCDIR="trilinos-$TRILINOS_VERSION-Source"

download_trilinos() {
    if [ -z "$TMPDIR" ]; then
        TMPDIR=/tmp
    fi
    TRILINOS_TARBALL=trilinos-$TRILINOS_VERSION-Source.tar.gz
    wget -P $TMPDIR -c http://trilinos.sandia.gov/download/files/$TRILINOS_TARBALL || exit 1
    tar -xzv -f $TMPDIR/$TRILINOS_TARBALL || exit 1
}

build_trilinos() {
    if [ ! -d "$TRILINOS_SRCDIR" ]; then
        download_trilinos
    fi
    cd $TRILINOS_SRCDIR || exit 1
    mkdir -p build
    cd build
    if [ -e CMakeCache.txt ]; then
        echo "Removing CMakeCache.txt"
        rm -f CMakeCache.txt
    fi
    $CMAKE \
        -Wno-dev \
        -D CMAKE_INSTALL_PREFIX:PATH="$PREFIX" \
        -D CMAKE_BUILD_TYPE:STRING=RELEASE \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -D TPL_ENABLE_MPI:BOOL=ON \
        -D TPL_ENABLE_Boost:BOOL=ON \
        -D TPL_BLAS_LIBRARIES="$BLAS_LIBS" \
        -D TPL_LAPACK_LIBRARIES="-lm $LAPACK_LIBS $BLAS_LIBS -lm" \
        -D Trilinos_ENABLE_OpenMP:BOOL=ON \
        -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=FALSE \
        -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=FALSE \
        -D Trilinos_ENABLE_Anasazi:BOOL=ON \
        -D Trilinos_ENABLE_Teuchos:BOOL=ON \
        $TRILINOS_EXTRA_CMAKE_OPTIONS \
        ../ || exit 1
    make || exit 1
    make install || exit 1
}

clean_trilinos() {
    target=clean
    if [ -n "$1" ]; then
        target="$1"
    fi
    cd $TRILINOS_SRCDIR || exit 1
    if [ -d build ]; then
        case $target in
        clean)
            cd build
            if [ -e Makefile ]; then
                make clean || exit 1
            fi
            ;;
        distclean)
            # Emulate make distclean
            rm -rf build/
            ;;
        uninstall)
            # Emulate make uninstall
            cd build
            if [ -e install_manifest.txt ]; then
                xargs rm -f < install_manifest.txt
            fi
            ;;
        esac
    fi
}
