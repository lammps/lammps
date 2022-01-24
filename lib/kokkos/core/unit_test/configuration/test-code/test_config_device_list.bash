
SRC_DIR=${KOKKOS_PATH}/core/unit_test/configuration/test-code
# List of parallel device types 
HostPDevices=(OpenMP Threads)
if [ ! -z "$KOKKOS_ARCH_TEST" ]; then
  HostPDevices=(OpenMP)
fi

if [ ! -z "$HPX_ROOT" ]
then 
  HostPDevices=(${HostPDevices[@]} HPX)
fi

if [ ! -z "$CUDA_ROOT" ]
then
  AccDevices=(${AccDevices[@]} Cuda)
  export CXX=${KOKKOS_PATH}/bin/nvcc_wrapper
fi
if [ ! -z "$HIP_ROOT" ]
then
  AccDevices=(${AccDevices[@]} HIP)
fi

for hpdevice in "${HostPDevices[@]}"
do
  hpdevice_up=`echo $hpdevice | tr a-z A-Z`
  CMAKE_HPDEVICE="-DKokkos_ENABLE_${hpdevice_up}=ON"

  if [ ! -z "$AccDevices" ]
  then
    for accdevice in "${AccDevices[@]}"
    do
      accdevice_up=`echo $accdevice | tr a-z A-Z`
      CMAKE_ACCDEVICE="-DKokkos_ENABLE_${accdevice_up}=ON"
      ${SRC_DIR}/test_config_arch_list.bash "$hpdevice,$accdevice" "${CMAKE_HPDEVICE} ${CMAKE_ACCDEVICE}"
      ${SRC_DIR}/test_config_arch_list.bash "$hpdevice,$accdevice,Serial" "${CMAKE_HPDEVICE} ${CMAKE_ACCDEVICE} -DKokkos_ENABLE_SERIAL=ON"
    done
  else
    #no, I need to be able to specify this
    #export CXX=g++
    ${SRC_DIR}/test_config_arch_list.bash "$hpdevice" "${CMAKE_HPDEVICE}"
    ${SRC_DIR}/test_config_arch_list.bash "$hpdevice,Serial" "${CMAKE_HPDEVICE} -DKokkos_ENABLE_SERIAL=ON"
  fi 
done

