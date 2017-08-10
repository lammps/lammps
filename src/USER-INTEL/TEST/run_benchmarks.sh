#!/bin/bash

#########################################################################
# Adjust settings below for your system
#########################################################################

# --------------------- MPI Launch Command

export MPI="mpirun"           
#export MPI="numactl -p 1 mpirun"    # -- Systems w/ MCDRAM in flat mode

# ------------- Name and location of the LAMMPS binary

export LMP_BIN=../../lmp_intel_cpu_intelmpi
#export LMP_BIN=../../lmp_knl

# ------------- Directory containing the LAMMPS installation

export LMP_ROOT=../../../

# ------------- Number of physical cores (not HW threads)

export LMP_CORES=36            # -- For Intel Xeon E5-2697v4 SKU
#export LMP_CORES=68           # -- For Intel Xeon Phi x200 7250 SKU

# ------------- Number of HW threads to use in tests

export LMP_THREAD_LIST="2"     # -- For 2 threads per core w/ HT enabled
#export LMP_THREAD_LIST="2 4"   # -- For 2 threads per core w/ HT enabled

# ------------- MPI Tuning Parameters

#export I_MPI_SHM_LMT=shm      # -- Uncomment for Xeon Phi x200 series

# ------------- Library locations for build

#source /opt/intel/parallel_studio_xe_2017.2.050/psxevars.sh

#########################################################################
# End settings for your system
#########################################################################

export WORKLOADS="lj rhodo rhodo_lrt lc sw water eam"
export LMP_ARGS="-pk intel 0 -sf intel -screen none -v d 1"
export RLMP_ARGS="-pk intel 0 lrt yes -sf intel -screen none -v d 1"

export LOG_DIR_HEADER=`echo $LMP_BIN | sed 's/\.\.\///g' | sed 's/\.\///g'`
export LOG_DIR_HOST=`hostname`
export DATE_STRING=`date +%s`
export LOG_DIR=$LOG_DIR_HOST"_"$LOG_DIR_HEADER"_"$DATE_STRING
mkdir $LOG_DIR

export I_MPI_PIN_DOMAIN=core
export I_MPI_FABRICS=shm
export KMP_BLOCKTIME=0

echo -n "Creating restart file...."
$MPI -np $LMP_CORES $LMP_BIN -in in.lc_generate_restart -log none $LMP_ARGS
echo "Done."
for threads in $LMP_THREAD_LIST
do
  export OMP_NUM_THREADS=$threads
  for workload in $WORKLOADS
  do
    export LOGFILE=$LOG_DIR/$workload.$LMP_CORES"c"$threads"t".log
    echo "Running $LOGFILE"
    cmd="$MPI -np $LMP_CORES $LMP_BIN -in in.intel.$workload -log $LOGFILE $LMP_ARGS";
    rthreads=$threads
    unset KMP_AFFINITY
    $cmd

    # - For benchmarks with PPPM, also try LRT mode
    if [ $workload = "rhodo" ]; then
      export LOGFILE=$LOG_DIR/$workload"_lrt".$LMP_CORES"c"$threads"t".log
      cmd="$MPI -np $LMP_CORES $LMP_BIN -in in.intel.$workload -log $LOGFILE $RLMP_ARGS";
      rthreads=$((threads-1))
      export KMP_AFFINITY=none
      export OMP_NUM_THREADS=$rthreads
      echo "  $cmd" >> $LOG_DIR/commands.info
      $cmd
    fi
  done
done

# Performance reported by LAMMPS (Timesteps/second ignoring warm-up run)
grep Perf $LOG_DIR/*.log | awk 'BEGIN{n=1}n%2==0{print $0}{n++}' | sed 's/\/day//g' | sed 's/steps\/s/steps_s/g' | sed 's/hours\/ns//g' | sed 's/.*\///g' | sed 's/\.log:Performance://g' | awk '{c=NF-1; print $1,$c}'
