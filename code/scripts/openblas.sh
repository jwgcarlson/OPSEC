# Shell functions for downloading/building/cleaning OpenBLAS source code.
#
# We assume these methods are being called from build.sh or clean.sh, so that
# setup.sh has been sourced and all relevant configuration variables defined.

# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    exit 1
fi

# Choose OpenBLAS version (default to latest git)
if [ -z "$OPENBLAS_VERSION" ]; then
    OPENBLAS_VERSION=git
fi
if [ -z "$OPENBLAS_SRCDIR" ]; then
    OPENBLAS_SRCDIR=OpenBLAS-$OPENBLAS_VERSION
fi
if [ -z "$OPENBLAS_TARBALL_URL" ]; then
    OPENBLAS_TARBALL_URL=http://github.com/xianyi/OpenBLAS/tarball
fi
if [ -z "$OPENBLAS_GIT_URL" ]; then
    OPENBLAS_GIT_URL=git://github.com/xianyi/OpenBLAS
fi


download_openblas() {
    if [ "$OPENBLAS_VERSION" = git ]; then
        GIT=${GIT:-git}
        if [ -z "`which $GIT`" ]; then
            echo "Cannot find git executable.  Please install git and try again."
            echo "If installed to a non-standard path, set the GIT environment variable."
            exit 1
        fi
        $GIT clone $OPENBLAS_GIT_URL $OPENBLAS_SRCDIR
    else
        if [ -z "$TMPDIR" ]; then
            TMPDIR=/tmp
        fi
        wget -P $TMPDIR -c $OPENBLAS_TARBALL_URL/$OPENBLAS_VERSION || exit 1
        tar -xzv -f $TMPDIR/xianyi-OpenBLAS-$OPENBLAS_VERSION*.tar.gz || exit 1
        mv xianyi-OpenBLAS* $OPENBLAS_SRCDIR
    fi
}

build_openblas() {
    if [ ! -d "$OPENBLAS_SRCDIR" ]; then
        download_openblas
    fi
    cd $OPENBLAS_SRCDIR || exit 1
    make || exit 1
    make install PREFIX=$PREFIX || exit 1
}

clean_openblas() {
    cd $OPENBLAS_SRCDIR || exit 1
    case "$1" in
    clean)
        make $target || exit 1
        ;;
    uninstall)
        rm -f $OPSEC_ROOT/lib/libopenblas*
        rm -f $OPSEC_ROOT/include/cblas.h
        rm -f $OPSEC_ROOT/include/f77blas.h
        rm -f $OPSEC_ROOT/include/openblas_config.h
        ;;
    *)
        echo "Doing nothing in $OPENBLAS_SRCDIR."
        ;;
    esac
}
