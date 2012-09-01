#!/bin/bash
#PBS -q euclid_batch
#PBS -l nodes=16:ppn=8:euclid
#PBS -l walltime=20:00:00
#PBS -V
#PBS -j oe
#PBS -o stdout
#PBS -N OPSEC

TOP=${OPSEC_ROOT}               # OPSEC install directory
NODES=1                         # Number of computational nodes available
CORES=4                         # Number of cores per node
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
        $1 -c $CONFIG 2>&1  | tee -a out.$1
    elif [ $# = 2 ]; then
        mpirun -n $2 $1 -c $CONFIG 2>&1  | tee -a out.$1
    elif [ $# = 3 ]; then
        mpirun -n $2 -npernode $3 $1 -c $CONFIG 2>&1  | tee -a out.$1
    else
        echo "Too many arguments to run."
    fi
    retval=${PIPESTATUS[0]}     # Exit code of cmd1 in a series of piped commands "cmd1 | cmd2 | ..."
    echo "$1 finished at `date`" | tee -a out.$1
    return $retval
}

if [ ! -e "gmock-gals.abn" ]; then
    run gmock || die
fi

#if [ ! -e "gals.abn" ]; then
#    $BIN/realmock input_galfile=gmock-gals.abn input_coordsys=cartesian L=3000 \
#        output_galfile=gals.abn output_coordsys=spherical output_ngals=100000 \
#        maskfile=mask.hpx radialfile=radial_profile.dat \
#        RMin=1200. RMax=1500. origin="1500 1500 1500" || die
#fi
#
#if [ ! -e "nbar.dat" ]; then
#    $BIN/makeradial gals.abn mask.hpx nbar.dat || die
#fi

run basis || die                # single node, OpenMP parallelism
run klt-anasazi || die         # hybrid OpenMP+MPI
run comma || die       # hybrid OpenMP+MPI
run dot || die                  # single node, OpenMP parallelism
run estimate || die
