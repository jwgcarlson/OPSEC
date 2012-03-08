#!/bin/bash
#
# Clean script for OPSEC.


die() { echo "Clean failed." && exit 1; }

# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    die
fi

# Determine absolute path to top-level source directory (the location of this script)
CODE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
    echo "Usage: $0 [-d | -u] [targets...]"
    echo "Valid clean targets: cfitsio, openblas, scalapack, arpack, trilinos, gmock, opsec"
    echo "Options:"
    echo "  -d     Run 'make distclean' instead of 'make clean'"
    echo "  -u     Run 'make uninstall' instead of 'make clean'"
}

# Load custom compilation flags
if [ -e "$CODE/setup.sh" ]; then
    source "$CODE/setup.sh"
fi

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

CLEANTARGET=clean

clean() {
    while [ -n "$1" ]; do
        cd $CODE
        case $1 in
            cfitsio)
                source $CODE/scripts/cfitsio.sh
                clean_cfitsio $CLEANTARGET || die
                ;;
            openblas)
                source $CODE/scripts/openblas.sh
                clean_openblas $CLEANTARGET || die
                ;;
            scalapack)
                source $CODE/scripts/scalapack.sh
                clean_scalapack $CLEANTARGET || die
                ;;
            arpack)
                cd $CODE/arpack
                echo "Cleaning all in $PWD"
                if [ -e Makefile ]; then
                    make $CLEANTARGET || die
                fi
                ;;
            trilinos)
                source $CODE/scripts/trilinos.sh
                clean_trilinos $CLEANTARGET || die
                ;;
            gmock)
                cd $CODE/gmock
                echo "Cleaning all in $PWD"
                if [ -e Makefile ]; then
                    make $CLEANTARGET || die
                fi
                ;;
            opsec)
                cd $CODE/opsec
                echo "Cleaning all in $PWD"
                if [ -e Makefile ]; then
                    make $CLEANTARGET || die
                fi
                ;;
            all)
                clean cfitsio arpack trilinos gmock opsec
                ;;
            *)
                echo "clean.sh: invalid clean target: $1"
                exit 1
                ;;
        esac
        shift
    done
}

while [ -n "$1" ]; do
    case $1 in
        -d)
            CLEANTARGET=distclean
            ;;
        -u)
            CLEANTARGET=uninstall
            ;;
        -*)
            echo "clean.sh: invalid command line switch: $1"
            exit 1
            ;;
        *)
            clean $1
            ;;
    esac
    shift
done

echo "Clean complete."
