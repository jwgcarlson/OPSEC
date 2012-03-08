# Shell functions for downloading/building/cleaning ScaLAPACK source code.
#
# We assume these methods are being called from build.sh or clean.sh, so that
# setup.sh has been sourced and all relevant configuration variables defined.

# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    exit 1
fi

# We need Python to run the ScaLAPACK installer
if [ -z "$OPSEC_PYTHON" ]; then
    echo "Python is required to install ScaLAPACK.  Please make sure a Python"
    echo "interpreter is installed, and set the OPSEC_PYTHON variable in opsecrc.sh."
    exit 1
fi

if [ -z "$SCALAPACK_INSTALLER_URL" ]; then
    SCALAPACK_INSTALLER_URL=http://www.netlib.org/scalapack
fi
if [ -z "$SCALAPACK_INSTALLER_TARBALL" ]; then
    SCALAPACK_INSTALLER_TARBALL=scalapack_installer.tgz
fi
if [ -z "$SCALAPACK_INSTALLER_DIR" ]; then
    SCALAPACK_INSTALLER_DIR=scalapack_installer
fi


download_scalapack_installer() {
    if [ -z "$TMPDIR" ]; then
        TMPDIR=/tmp
    fi
    wget -P $TMPDIR -c $SCALAPACK_INSTALLER_URL/$SCALAPACK_INSTALLER_TARBALL || exit 1
    tar -xzv -f $TMPDIR/$SCALAPACK_INSTALLER_TARBALL || exit 1
    mv scalapack_installer_* $SCALAPACK_INSTALLER_DIR || exit 1
}

build_scalapack() {
    if [ ! -d "$SCALAPACK_INSTALLER_DIR" ]; then
        download_scalapack_installer
    fi
    cd $SCALAPACK_INSTALLER_DIR || exit 1
    options=(--prefix="$PREFIX" --notesting --mpiincdir=none)
    if [ -n "$MPICC" ]; then
        options+=(--mpicc=$MPICC)
    fi
    if [ -n "$MPIF77" ]; then
        options+=(--mpicc=$MPIF77)
    fi
    if [ -n "$CFLAGS" ]; then
        options+=(--ccflags="$CFLAGS")
    fi
    if [ -n "$FFLAGS" ]; then
        options+=(--fcflags="$FFLAGS")
    fi
    if [ -n "$BLAS_LIBS" ]; then
        options+=(--blaslib="$LDFLAGS $BLAS_LIBS")
    fi
    $OPSEC_PYTHON setup.py "${options[@]}"
}

clean_scalapack() {
    cd $SCALAPACK_INSTALLER_DIR || exit 1
    case "$1" in
    clean|distclean)
        $OPSEC_PYTHON setup.py --clean || exit 1
        ;;
    uninstall)
        rm -f $PREFIX/lib/libscalapack.a
        ;;
    esac
}
