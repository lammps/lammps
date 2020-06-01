SRC_DIR=${KOKKOS_PATH}/core/unit_test/configuration/test-code

# List of parallel device types 
Options=(deprecated_code aggressive_vectorization disable_profiling large_mem_tests)
CudaOptions=(lambda relocatable_device_code uvm constexpr)

if [ ! -z "$KOKKOS_ARCH_TEST" ]; then
  Options=(disable_profiling)
  CudaOptions=(uvm)
fi

MakeDevices=$1
CMakeDevices=$2
MakeArch=$3
CMakeArch=$4

for option in "${Options[@]}"
do
  option_up=`echo $option | tr a-z A-Z`
  if [[ $option_up == *"DISABLE"* ]]; then
    new_option_up=${option_up/DISABLE_/}
    CMAKE_OPTION="-DKokkos_ENABLE_${new_option_up}=OFF"
  else
    CMAKE_OPTION="-DKokkos_ENABLE_${option_up}=ON"
  fi

  #Renaming options as GNU Make expects them
  option=${option/deprecated_code/enable_deprecated_code}
  option=${option/large_mem_tests/enable_large_mem_tests}

  if [ ! -z $CudaOptions ]; then 
    for cuda_option in "${CudaOptions[@]}"
    do
      cuda_option_up=`echo $cuda_option | tr a-z A-Z`
      CMAKE_CUDA_OPTION="-DKokkos_ENABLE_CUDA_${cuda_option_up}=ON"

      #Renaming options as GNU Make expects them
      cuda_option=${cuda_option/lambda/enable_lambda}
      cuda_option=${cuda_option/constexpr/enable_constexpr}
      cuda_option=${cuda_option/relocatable_device_code/rdc}
      cuda_option=${cuda_option/uvm/force_uvm}

      ${SRC_DIR}/test_config_run.bash "$MakeDevices" "$CMakeDevices" "$MakeArch" "$CMakeArch" "KOKKOS_OPTIONS=$option KOKKOS_CUDA_OPTIONS=$cuda_option" "$CMAKE_OPTION $CMAKE_CUDA_OPTION"
    done
  else  
    ${SRC_DIR}/test_config_run.bash "$MakeDevices" "$CMakeDevices" "$MakeArch" "$CMakeArch" "KOKKOS_OPTIONS=$option" "$CMAKE_OPTION"
  fi
done

