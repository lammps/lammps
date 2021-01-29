#!/bin/bash

# BytesAndFlops
cd build/bytes_and_flops

USE_CUDA=`grep "_CUDA" KokkosCore_config.h | wc -l`

if [[ ${USE_CUDA} > 0 ]]; then
  BAF_EXE=bytes_and_flops.cuda
  TEAM_SIZE=256
else
  BAF_EXE=bytes_and_flops.exe
  TEAM_SIZE=1
fi

BAF_PERF_1=`./${BAF_EXE} 2 100000 1024 1 1 1 1 ${TEAM_SIZE} 6000 | awk '{print $12/174.5}'`
BAF_PERF_2=`./${BAF_EXE} 2 100000 1024 16 1 8 64 ${TEAM_SIZE} 6000 | awk '{print $14/1142.65}'`

echo "BytesAndFlops: ${BAF_PERF_1} ${BAF_PERF_2}"
cd ../..


# MiniMD
cd build/miniMD
cp ../../miniMD/kokkos/Cu_u6.eam ./
MD_PERF_1=`./miniMD --half_neigh 0 -s 60 --ntypes 1 -t ${OMP_NUM_THREADS} -i ../../miniMD/kokkos/in.eam.miniMD | grep PERF_SUMMARY | awk '{print $10/21163341}'`
MD_PERF_2=`./miniMD --half_neigh 0 -s 20 --ntypes 1 -t ${OMP_NUM_THREADS} -i ../../miniMD/kokkos/in.eam.miniMD | grep PERF_SUMMARY | awk '{print $10/13393417}'`

echo "MiniMD: ${MD_PERF_1} ${MD_PERF_2}"
cd ../..

# MiniFE
cd build/miniFE
rm *.yaml
./miniFE.x -nx 100 &> /dev/null
FE_PERF_1=`grep "CG Mflop" *.yaml | awk '{print $4/14174}'`
rm *.yaml
./miniFE.x -nx 50 &> /dev/null
FE_PERF_2=`grep "CG Mflop" *.yaml | awk '{print $4/11897}'`
cd ../..
echo "MiniFE: ${FE_PERF_1} ${FE_PERF_2}"

PERF_RESULT=`echo "${BAF_PERF_1} ${BAF_PERF_2} ${MD_PERF_1} ${MD_PERF_2} ${FE_PERF_1} ${FE_PERF_2}" | awk '{print ($1+$2+$3+$4+$5+$6)/6}'`
echo "Total Result: " ${PERF_RESULT}
