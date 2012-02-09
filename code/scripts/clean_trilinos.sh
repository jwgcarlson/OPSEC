#!/bin/sh
#
# Script for cleaning Trilinos source tree.
#
# We assume that this is being called from OPSEC's clean.sh script, so that
# all necessary environment variables are already defined.


# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    die
fi

die() {
    echo "Trilinos clean failed."
    exit 1
}

# Choose Trilinos version
if [ -z "$TRILINOS_VERSION" ]; then
    TRILINOS_VERSION=10.8.5
fi
TRILINOS_SRCDIR="trilinos-$TRILINOS_VERSION-Source"

# Clean
echo "$PWD"
cd $TRILINOS_SRCDIR || die
echo "Cleaning all in $PWD"
if [ -d build ]; then
    cd build
    if [ -e Makefile ]; then
        make clean || die
    fi
fi
