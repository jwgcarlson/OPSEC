#!/bin/bash
#PBS -q euclid_batch
#PBS -l nodes=16:ppn=8:euclid
#PBS -l walltime=16:00:00
#PBS -V
#PBS -j oe
#PBS -o stdout
#PBS -N OPSEC

TOP=${TOP:-$HOME/gpower}        # Top-level gpower directory
BIN=$TOP/bin                    # gpower/bin directory
NODES=16                        # Number of computational nodes available
CORES=8                         # Number of cores per node
PROCESSES=$(($NODES*$CORES))    # Number of MPI processes for fully parallel job
CONFIG=opsec.cfg                # Config file

# Change to work directory if this is a batch job
if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi

die() {
    echo "Run failed."
    exit 1
}

run() {
    echo "Starting $1 at `date`" | tee out.$1
    if [ $# = 1 ]; then
        $BIN/$1 -c $CONFIG 2>&1  | tee -a out.$1
    elif [ $# = 2 ]; then
        mpirun -n $2 $BIN/$1 -c $CONFIG 2>&1  | tee -a out.$1
    elif [ $# = 3 ]; then
        mpirun -n $2 -npernode $3 $BIN/$1 -c $CONFIG 2>&1  | tee -a out.$1
    else
        echo "Too many arguments to run."
    fi
    retval=${PIPESTATUS[0]}     # Exit code of cmd1 in a series of piped commands "cmd1 | cmd2 | ..."
    echo "$1 finished at `date`" | tee -a out.$1
    return $retval
}

export OMP_NUM_THREADS=$CORES

run basis || die                # single node, OpenMP parallelism
run klt $NODES 1 || die         # hybrid OpenMP+MPI
run comma $NODES 1 || die       # hybrid OpenMP+MPI
run dot || die                  # single node, OpenMP parallelism
run estimate || die
