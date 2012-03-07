# Shell functions for downloading/building/cleaning Cuba source code.
#
# We assume these methods are being called from build.sh or clean.sh, so that
# setup.sh has been sourced and all relevant configuration variables defined.

# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    exit 1
fi

# Choose Cuba version (default to 3.0)
if [ -z "$CUBA_VERSION" ]; then
    CUBA_VERSION=3.0
fi

CUBA_SRCDIR="Cuba-$CUBA_VERSION"

download_cuba() {
    if [ -z "$TMPDIR" ]; then
        TMPDIR=/tmp
    fi
    CUBA_TARBALL=Cuba-$CUBA_VERSION.tar.gz
    wget -P $TMPDIR -c http://www.feynarts.de/cuba/$CUBA_TARBALL || exit 1
    tar -xzv -f $TMPDIR/$CUBA_TARBALL || exit 1
}

build_cuba() {
    if [ ! -d "$CUBA_SRCDIR" ]; then
        download_cuba
    fi
    cd $CUBA_SRCDIR || exit 1
    ./configure --prefix="$PREFIX" || exit 1
    make || exit 1
    make install || exit 1
}

clean_cuba() {
    target=clean
    if [ -n "$1" ]; then
        target="$1"
    fi
    cd $CUBA_SRCDIR || exit 1
    if [ -e makefile ]; then
        make $target || exit 1
    fi
}
