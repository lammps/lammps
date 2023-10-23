#!/bin/bash

#########################################################################
# Adjust settings below for your system
#########################################################################

# --------------------- MPI Launch Command

export MPI="mpirun"

# ------------- Name and location of the LAMMPS binary

export LMP_BIN=../../lmp_intel_cpu_intelmpi

# ------------- Directory containing the LAMMPS installation

export LMP_ROOT=../../../

# ------------- Number of physical cores (not HW threads)

export LMP_CORES=36
# Set automatically with lscpu
export LMP_CORES=`lscpu | awk '$1=="Core(s)"{t=NF; cores=$t}$1=="Socket(s):"{t=NF; sockets=$t}END{print cores*sockets}'`

# ------------- Number of HW threads to use in tests

export LMP_THREAD_LIST="1 2"   # -- Also test 2 threads per core w/ HT enabled

# ------------- MPI Tuning Parameters

export I_MPI_FABRICS=shm
export I_MPI_PIN_DOMAIN=core

#########################################################################
# End settings for your system
#########################################################################

export WORKLOADS="lj rhodo lc sw water eam airebo dpd tersoff snap"
export LMP_ARGS="-pk intel 0 -sf intel -screen none -v d 1"
export RLMP_ARGS="-pk intel 0 lrt yes -sf intel -screen none -v d 1"

export LOG_DIR_HEADER=`echo $LMP_BIN | sed 's/\.\.\///g' | sed 's/\.\///g'`
export LOG_DIR_HOST=`hostname`
export DATE_STRING=`date +%s`
export LOG_DIR=$LOG_DIR_HOST"_"$LOG_DIR_HEADER"_"$DATE_STRING
mkdir $LOG_DIR

export KMP_BLOCKTIME=0

echo -n "Creating restart file...."
$MPI -np $LMP_CORES $LMP_BIN -in in.lc_generate_restart -log none $LMP_ARGS
echo "Done."
for threads in $LMP_THREAD_LIST
do
  for workload in $WORKLOADS
  do
    export LOGFILE=$LOG_DIR/$workload.$LMP_CORES"c"$threads"t".log
    echo "Running $LOGFILE"
    cmd="$MPI -np $LMP_CORES $LMP_BIN -in in.intel.$workload -log $LOGFILE $LMP_ARGS";
    export OMP_NUM_THREADS=$threads
    unset KMP_AFFINITY
    $cmd

    # - For benchmarks with PPPM, also try LRT mode
    if [ $workload = "rhodo" ]; then
      export LOGFILE=$LOG_DIR/$workload"_lrt".$LMP_CORES"c"$threads"t".log
      cmd="$MPI -np $LMP_CORES $LMP_BIN -in in.intel.$workload -log $LOGFILE $RLMP_ARGS";
      rthreads=$((threads-1))
      if [ $rthreads -lt 1 ]; then
          rthreads=1
      fi
      export KMP_AFFINITY=none
      export OMP_NUM_THREADS=$rthreads
      echo "  $cmd" >> $LOG_DIR/commands.info
      $cmd
    fi
  done
done

# Performance reported by LAMMPS (Timesteps/second ignoring warm-up run)
grep Perf $LOG_DIR/*.log | awk 'n%2==1{c=NF-3; print $1,$c}{n++}' | sed -e 's/.*\///g' -e 's/:.*://g'

