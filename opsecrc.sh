# OPSEC environment definitions for bash (and related) shells.  Modify as
# appropriate.
#
# Source this file before each OPSEC session.  Alternatively, source it in your
# .bashrc file.


# Location of this file.
thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Top-level OPSEC install directory.  This default can be overridden.
OPSEC_ROOT="$thisdir"

# Choose a Python interpreter.  Leave blank to disable Python support.
OPSEC_PYTHON=python

# Verify Python interpreter is present and working.  (Python 2 only!)
if [ -n "$OPSEC_PYTHON" ]; then
    if [ -z "`which $OPSEC_PYTHON`" -o "`$OPSEC_PYTHON -c 'print "x"' 2>/dev/null`" != x ]; then
        echo "** Python interpreter '$OPSEC_PYTHON' not found.  Python support disabled."
        OPSEC_PYTHON=""
    else
        PYTHON_VERSION=`$OPSEC_PYTHON -c "import sys; print sys.version[:3]"`
    fi
fi


# Export variables and update paths.
export OPSEC_ROOT
export PATH="$PATH:$OPSEC_ROOT/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$OPSEC_ROOT/lib"
if [ -n "$OPSEC_PYTHON" ]; then
    export OPSEC_PYTHON
    export PYTHONPATH="$PYTHONPATH:$OPSEC_ROOT/lib/python$PYTHON_VERSION"
fi
