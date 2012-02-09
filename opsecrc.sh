# OPSEC environment definitions for bash (and related) shells.  Modify as
# appropriate.
#
# Source this file before each OPSEC session.  Alternatively, source it in your
# .bashrc file.


# Location of this file.
cd `dirname $0`
THISDIR="$PWD"

# Top-level OPSEC install directory.  This default can be overridden.
OPSEC_ROOT="$THISDIR"

# Choose a Python interpreter.  Leave blank to disable Python support.
OPSEC_PYTHON=python

# Verify Python interpreter is present and working.  (Python 2 only!)
if [ -n "$OPSEC_PYTHON" ]; then
    if [ -z "`which $OPSEC_PYTHON`" -o "`$OPSEC_PYTHON -c 'print "x"' 2>/dev/null`" != x ]; then
        echo "** Python interpreter '$OPSEC_PYTHON' not found.  Python support disabled."
        OPSEC_PYTHON=""
    fi
fi

# Export variables and update paths.
export OPSEC_ROOT
export PATH="$PATH:$OPSEC_ROOT/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$OPSEC_ROOT/lib"
if [ -n "$OPSEC_PYTHON" ]; then
    export OPSEC_PYTHON
    export PYTHONPATH="$PYTHONPATH:$OPSEC_ROOT/lib/python"
fi
