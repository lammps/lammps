if(CMAKE_PROJECT_NAME STREQUAL "Trilinos")
  set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Whether to build Serial backend" FORCE)

  if(NOT ${Trilinos_ENABLE_OpenMP} STREQUAL "")
    set(Kokkos_ENABLE_OPENMP ${Trilinos_ENABLE_OpenMP} CACHE BOOL "Whether to build OpenMP backend" FORCE)
  else()
    set(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "Whether to build OpenMP backend" FORCE)
  endif()

  if(NOT ${TPL_ENABLE_CUDA} STREQUAL "")
    set(Kokkos_ENABLE_CUDA ${TPL_ENABLE_CUDA} CACHE BOOL "Whether to build CUDA backend" FORCE)
  else()
    set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "Whether to build CUDA backend" FORCE)
  endif()

  if(NOT ${TPL_ENABLE_HPX} STREQUAL "")
    set(Kokkos_ENABLE_HPX ${TPL_ENABLE_HPX} CACHE BOOL "Whether to build HPX backend" FORCE)
  else()
    set(Kokkos_ENABLE_HPX OFF CACHE BOOL "Whether to build HPX backend" FORCE)
  endif()

  if(NOT ${TPL_ENABLE_quadmath} STREQUAL "")
    set(Kokkos_ENABLE_LIBQUADMATH ${TPL_ENABLE_quadmath} CACHE BOOL "Whether to enable the LIBQUADMATH library" FORCE)
  else()
    set(Kokkos_ENABLE_LIBQUADMATH OFF CACHE BOOL "Whether to enable the LIBQUADMATH library" FORCE)
  endif()

  if(NOT ${TPL_ENABLE_DLlib} STREQUAL "")
    set(Kokkos_ENABLE_LIBDL ${TPL_ENABLE_DLlib} CACHE BOOL "Whether to enable the LIBDL library" FORCE)
  else()
    set(Kokkos_ENABLE_LIBDL OFF CACHE BOOL "Whether to enable the LIBDL library" FORCE)
  endif()

  set(Kokkos_ENABLE_COMPLEX_ALIGN OFF CACHE BOOL "Whether to align Kokkos::complex to 2*alignof(RealType)")

  # FIXME_TRILINOS We run into problems when trying to use an external GTest in Trilinos CI
  set(CMAKE_DISABLE_FIND_PACKAGE_GTest ON)

  assert_defined(
    Kokkos_ENABLE_SERIAL
    Kokkos_ENABLE_OPENMP
    Kokkos_ENABLE_THREADS
    Kokkos_ENABLE_HPX
    Kokkos_ENABLE_CUDA
    Kokkos_ENABLE_HIP
    Kokkos_ENABLE_SYCL
    Kokkos_ENABLE_OPENACC
  )

  if(Kokkos_ENABLE_SERIAL
     AND NOT Kokkos_ENABLE_OPENMP
     AND NOT Kokkos_ENABLE_THREADS
     AND NOT Kokkos_ENABLE_HPX
     AND NOT Kokkos_ENABLE_CUDA
     AND NOT Kokkos_ENABLE_HIP
     AND NOT Kokkos_ENABLE_SYCL
     AND NOT Kokkos_ENABLE_OPENACC
  )
    if(NOT DEFINED Kokkos_ENABLE_ATOMICS_BYPASS)
      set(Kokkos_ENABLE_ATOMICS_BYPASS ON
          CACHE BOOL "Whether to make atomics non-atomic for non-threaded MPI-only use cases" FORCE
      )
    endif()
  endif()
endif()
