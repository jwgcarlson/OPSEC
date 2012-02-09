#!/bin/bash
#
# Build script for OPSEC.


# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    die
fi

# Determine absolute path to top-level source directory (the location of this script)
cd `dirname $0`
export TOP=$PWD

die() {
    echo "Build failed."
    exit 1
}

usage() {
    echo "Usage: $0 all"
    echo "       $0 [targets...]"
    echo "Valid build targets: cfitsio, arpack, trilinos, gmock, opsec"
    echo ""
    echo "Note that dependencies are not respected when building targets individually."
    echo "You must make sure to build targets in the correct order, in accordance with"
    echo "the dependencies below.  (External or optional dependencies in parentheses.)"
    echo ""
    echo "Library dependencies:"
    echo "  cfitsio: none"
    echo "  arpack: (BLAS), (LAPACK)"
    echo "  trilinos: (BLAS), (LAPACK), (Boost)"
    echo "  gmock: (FFTW 2)"
    echo "  opsec: cfitsio, arpack, trilinos, (BLAS), (CImg.h), (libpng)"
}

# Load custom compilation flags
if [ -e "$TOP/setup.sh" ]; then
    source "$TOP/setup.sh"
fi

# Default installation prefix
if [ -z "$PREFIX" ]; then
    export PREFIX="$OPSEC_ROOT"
fi

# Set initial compiler and linker search directories
export CPPFLAGS="-I$PREFIX/include $CPPFLAGS"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS"

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

build() {
    while [ -n "$1" ]; do
        case $1 in
            cfitsio)
                cd $TOP
                $TOP/scripts/build_cfitsio.sh || die
                ;;
            arpack)
                cd $TOP/arpack
                echo "Building all in $PWD"
                ./configure --prefix=$PREFIX && make && make install || die
                ;;
            trilinos)
                cd $TOP
                $TOP/scripts/build_trilinos.sh || die
                ;;
            gmock)
                cd $TOP/gmock
                echo "Building all in $PWD"
                ./configure --prefix=$PREFIX && make && make install || die
                ;;
            opsec)
                cd $TOP/opsec
                echo "Building all in $PWD"
                ./configure --prefix=$PREFIX && make && make install || die
                ;;
            all)
                build cfitsio arpack trilinos gmock opsec
                ;;
            *)
                echo "build.sh: invalid build target: $1"
                exit 1
                ;;
        esac
        shift
    done
}

while [ -n "$1" ]; do
    case $1 in
        -h)
            usage
            ;;
        -*)
            echo "build.sh: invalid command line switch: $1"
            exit 1
            ;;
        *)
            build $1
            ;;
    esac
    shift
done

echo "Build complete."
