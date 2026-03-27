function(KOKKOS_DEVICE_OPTION SUFFIX DEFAULT DEV_TYPE DOCSTRING)
  kokkos_option(ENABLE_${SUFFIX} ${DEFAULT} BOOL ${DOCSTRING})
  string(TOUPPER ${SUFFIX} UC_NAME)
  if(KOKKOS_ENABLE_${UC_NAME})
    list(APPEND KOKKOS_ENABLED_DEVICES ${SUFFIX})
    #I hate that CMake makes me do this
    set(KOKKOS_ENABLED_DEVICES ${KOKKOS_ENABLED_DEVICES} PARENT_SCOPE)
  endif()
  set(KOKKOS_ENABLE_${UC_NAME} ${KOKKOS_ENABLE_${UC_NAME}} PARENT_SCOPE)
  if(KOKKOS_ENABLE_${UC_NAME} AND DEV_TYPE STREQUAL "HOST")
    set(KOKKOS_HAS_HOST ON PARENT_SCOPE)
  endif()
endfunction()

kokkos_cfg_depends(DEVICES NONE)

# Put a check in just in case people are using this option
kokkos_deprecated_list(DEVICES ENABLE)

kokkos_device_option(THREADS OFF HOST "Whether to build C++ threads backend")

# detect clang++ / cl / clang-cl clashes
if(CMAKE_CXX_COMPILER_ID STREQUAL Clang AND "x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
  # this specific test requires CMake >= 3.15
  if("x${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}" STREQUAL "xGNU")
    # use pure clang++ instead of clang-cl
    set(KOKKOS_COMPILER_CLANG_MSVC OFF)
  else()
    # it defaults to clang-cl
    set(KOKKOS_COMPILER_CLANG_MSVC ON)
  endif()
endif()

if(Trilinos_ENABLE_Kokkos AND Trilinos_ENABLE_OpenMP)
  set(OMP_DEFAULT ON)
else()
  set(OMP_DEFAULT OFF)
endif()
kokkos_device_option(OPENMP ${OMP_DEFAULT} HOST "Whether to build OpenMP backend")

# We want this to default to OFF for cache reasons, but if no
# host space is given, then activate serial
if(KOKKOS_HAS_HOST)
  set(SERIAL_DEFAULT OFF)
else()
  set(SERIAL_DEFAULT ON)
  if(NOT DEFINED Kokkos_ENABLE_SERIAL)
    message(
      STATUS
        "SERIAL backend is being turned on to ensure there is at least one Host space. To change this, you must enable another host execution space and configure with -DKokkos_ENABLE_SERIAL=OFF."
    )
  endif()
endif()
kokkos_device_option(SERIAL ${SERIAL_DEFAULT} HOST "Whether to build serial backend")

kokkos_device_option(HPX OFF HOST "Whether to build HPX backend (experimental)")

# Device backends have to come after host backends for header include order reasons
# Without this we can't make e.g. CudaSpace accessible by HostSpace
kokkos_device_option(OPENACC OFF DEVICE "Whether to build the OpenACC backend")
if(KOKKOS_ENABLE_OPENACC)
  compiler_specific_flags(
    Clang
    -fopenacc
    -fopenacc-fake-async-wait
    -fopenacc-implicit-worker=vector
    -Wno-openacc-and-cxx
    -Wno-openmp-mapping
    -Wno-unknown-cuda-version
    -Wno-pass-failed
  )
  compiler_specific_defs(Clang KOKKOS_WORKAROUND_OPENMPTARGET_CLANG)
endif()

kokkos_device_option(OPENMPTARGET OFF DEVICE "Whether to build the OpenMP target backend")
if(KOKKOS_ENABLE_OPENMPTARGET)
  set(ClangOpenMPFlag -fopenmp=libomp)

  compiler_specific_flags(
    Clang
    ${ClangOpenMPFlag}
    -Wno-openmp-mapping
    -Wno-unknown-cuda-version
    -Wno-pass-failed
    DEFAULT
    -fopenmp
  )
  compiler_specific_defs(Clang KOKKOS_WORKAROUND_OPENMPTARGET_CLANG)
endif()

if(Trilinos_ENABLE_Kokkos AND TPL_ENABLE_CUDA)
  set(CUDA_DEFAULT ON)
else()
  set(CUDA_DEFAULT OFF)
endif()
kokkos_device_option(CUDA ${CUDA_DEFAULT} DEVICE "Whether to build CUDA backend")

if(KOKKOS_ENABLE_CUDA)
  global_set(KOKKOS_DONT_ALLOW_EXTENSIONS "CUDA enabled")
  ## Cuda has extra setup requirements, turn on Kokkos_Setup_Cuda.hpp in macros
  list(APPEND DEVICE_SETUP_LIST Cuda)
endif()

kokkos_device_option(HIP OFF DEVICE "Whether to build HIP backend")

## HIP has extra setup requirements, turn on Kokkos_Setup_HIP.hpp in macros
if(KOKKOS_ENABLE_HIP)
  list(APPEND DEVICE_SETUP_LIST HIP)
endif()

kokkos_device_option(SYCL OFF DEVICE "Whether to build SYCL backend")

## SYCL has extra setup requirements, turn on Kokkos_Setup_SYCL.hpp in macros
if(KOKKOS_ENABLE_SYCL)
  list(APPEND DEVICE_SETUP_LIST SYCL)
endif()
