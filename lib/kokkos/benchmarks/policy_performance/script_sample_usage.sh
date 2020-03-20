#!/bin/bash

# Sample script for benchmarking policy performance 

# Suggested environment variables to export prior to executing script:
# KNL: 
# OMP_NUM_THREADS=256 KMP_AFFINITY=compact
# Power:
# OMP_NUM_THREADS=64 OMP_PROC_BIND=true

# Constants and Variables:
# Vary:  TEAMSIZE, and THREADRANGE
#  for TEAMSIZE in {1,2,4,5,8}; do
#  for THREADRANGE in {32,41,1000}; do
# Fixed: TEAMRANGE, VECTORRANGE, VECTORSIZE
# System specific: Adjust REPEAT values to architecture tests are run on

# Tests
# Static SCHEDULE = 1
# Tier 1: parallel_for + RangePolicy 300
# Tier 2: parallel_reduce, parallel_scan + RangePolicy 400 500
# Tier 3: 'outer' parallel_for with TeamPolicy (nested parallelism) 1XY
# Tier 4: 'outer' parallel_reduce with TeamPolicy (nested parallelism) 2XY
# Dynamic SCHEDULE = 2
# Tier 5: parallel_for + RangePolicy 300
# Tier 6: parallel_reduce, parallel_scan + RangePolicy 400 500
# Tier 7: 'outer' parallel_for with TeamPolicy (nested parallelism) 1XY
# Tier 8: 'outer' parallel_reduce with TeamPolicy (nested parallelism) 2XY

# Results grouped by: 
# 0) SCHEDULE  1) CODE (test)  2) TEAMRANGE  3) TEAMSIZE  4) THREADRANGE

EXECUTABLE=policy_performance

# Default defined values
TEAMRANGE=1000
THREADRANGE=1
VECTORRANGE=32
TEAMSIZE=1
VECTORSIZE=1
OREPEAT=1
MREPEAT=1
IREPEAT=1
SCHEDULE=1

# Host tests
SUFFIX=host
if [ -e $EXECUTABLE.$SUFFIX ]; then
echo "Host"

for SCHEDULE in {1,2}; do

# Tier 1 and 2, 5 and 6
for CODE in {300,400,500}; do
    for TEAMSIZE in {1,2,4,5,8}; do
    OMP_PROC_BIND=true ./$EXECUTABLE.$SUFFIX $TEAMRANGE $THREADRANGE $VECTORRANGE $OREPEAT $MREPEAT $IREPEAT $TEAMSIZE $VECTORSIZE $SCHEDULE $CODE
    done
done

# Tier 3, 7
for CODE in {100,110,111,112,120,121,122}; do
    for TEAMSIZE in {1,2,4,5,8}; do
      for THREADRANGE in {32,41,1000}; do
      OMP_PROC_BIND=true ./$EXECUTABLE.$SUFFIX $TEAMRANGE $THREADRANGE $VECTORRANGE $OREPEAT $MREPEAT $IREPEAT $TEAMSIZE $VECTORSIZE $SCHEDULE $CODE
      done
    done
done

# Tier 4, 8
for CODE in {200,210,211,212,220,221,222}; do
    for TEAMSIZE in {1,2,4,5,8}; do
      for THREADRANGE in {32,41,1000}; do
      OMP_PROC_BIND=true ./$EXECUTABLE.$SUFFIX $TEAMRANGE $THREADRANGE $VECTORRANGE $OREPEAT $MREPEAT $IREPEAT $TEAMSIZE $VECTORSIZE $SCHEDULE $CODE
      done
    done
done

done # end SCHEDULE

fi # end host


# Cuda tests
SUFFIX=cuda
# TEAMRANGE=10000, TEAMSIZE=8 too large
# TEAMRANGE=10000, TEAMSIZE=8, THREADRANGE=1000 too large
if [ -e $EXECUTABLE.$SUFFIX ]; then
echo "Cuda"

for SCHEDULE in {1,2}; do

# Reset defaults
TEAMRANGE=1000
THREADRANGE=1
VECTORRANGE=32
TEAMSIZE=1
VECTORSIZE=1

# Tier 1 and 2, 5 and 6
for CODE in {300,400,500}; do
    for TEAMSIZE in {1,2,4,5,8}; do
    ./$EXECUTABLE.$SUFFIX $TEAMRANGE $THREADRANGE $VECTORRANGE $OREPEAT $MREPEAT $IREPEAT $TEAMSIZE $VECTORSIZE $SCHEDULE $CODE
    done
done

# Tier 3, 7
for CODE in {100,110,111,112,120,121,122}; do
    for TEAMSIZE in {1,2,4,5,8}; do
      for THREADRANGE in {32,41,1000}; do
      ./$EXECUTABLE.$SUFFIX $TEAMRANGE $THREADRANGE $VECTORRANGE $OREPEAT $MREPEAT $IREPEAT $TEAMSIZE $VECTORSIZE $SCHEDULE $CODE
      done
    done
done

# Tier 4, 8
for CODE in {200,210,211,212,220,221,222}; do
    for TEAMSIZE in {1,2,4,5,8}; do
      for THREADRANGE in {32,41,1000}; do
      ./$EXECUTABLE.$SUFFIX $TEAMRANGE $THREADRANGE $VECTORRANGE $OREPEAT $MREPEAT $IREPEAT $TEAMSIZE $VECTORSIZE $SCHEDULE $CODE
      done
    done
done

done # end SCHEDULE

fi #end cuda
