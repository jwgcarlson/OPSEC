#!/bin/bash
#
# Clean script for OPSEC.


# Ensure opsecrc.sh has been loaded
if [ -z "$OPSEC_ROOT" ]; then
    echo "You must first load the OPSEC environment by sourcing opsecrc.sh."
    die
fi

# Determine absolute path to top-level source directory (the location of this script)
cd `dirname $0`
export TOP=$PWD

die() {
    echo "Clean failed."
    exit 1
}

usage() {
    echo "Usage: $0 [-d] [targets...]"
    echo "Valid clean targets: cfitsio, arpack, trilinos, gmock, opsec"
    echo "Options:"
    echo "  -d     Run 'make distclean' instead of 'make clean'"
}

# Load custom compilation flags
if [ -e "$TOP/setup.sh" ]; then
    source "$TOP/setup.sh"
fi

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

CLEANTARGET=clean

clean() {
    while [ -n "$1" ]; do
        case $1 in
            cfitsio)
                cd $TOP/cfitsio
                echo "Cleaning all in $PWD"
                if [ -e Makefile ]; then
                    make clean || die
                fi
                ;;
            arpack)
                cd $TOP/arpack
                echo "Cleaning all in $PWD"
                if [ -e Makefile ]; then
                    make $CLEANTARGET || die
                fi
                ;;
            trilinos)
                cd $TOP
                $TOP/scripts/clean_trilinos.sh || die
                ;;
            gmock)
                cd $TOP/gmock
                echo "Cleaning all in $PWD"
                if [ -e Makefile ]; then
                    make $CLEANTARGET || die
                fi
                ;;
            opsec)
                cd $TOP/opsec
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
