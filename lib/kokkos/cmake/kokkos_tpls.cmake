kokkos_cfg_depends(TPLS OPTIONS)
kokkos_cfg_depends(TPLS DEVICES)
kokkos_cfg_depends(TPLS COMPILER_ID)

function(KOKKOS_TPL_OPTION PKG DEFAULT)
  cmake_parse_arguments(PARSED "" "TRIBITS" "" ${ARGN})

  if(PARSED_TRIBITS)
    #this is also a TPL option you can activate with Tribits
    if(NOT "${TPL_ENABLE_${PARSED_TRIBITS}}" STREQUAL "")
      #Tribits brought its own default that should take precedence
      set(DEFAULT ${TPL_ENABLE_${PARSED_TRIBITS}})
    endif()
  endif()

  kokkos_enable_option(${PKG} ${DEFAULT} "Whether to enable the ${PKG} library")
  kokkos_option(${PKG}_DIR "" PATH "Location of ${PKG} library")
  set(KOKKOS_ENABLE_${PKG} ${KOKKOS_ENABLE_${PKG}} PARENT_SCOPE)
  set(KOKKOS_${PKG}_DIR ${KOKKOS_${PKG}_DIR} PARENT_SCOPE)
endfunction()

kokkos_tpl_option(HWLOC Off TRIBITS HWLOC)
kokkos_tpl_option(CUDA ${Kokkos_ENABLE_CUDA} TRIBITS CUDA)
if(KOKKOS_ENABLE_HIP AND NOT KOKKOS_CXX_COMPILER_ID STREQUAL HIPCC)
  set(ROCM_DEFAULT ON)
else()
  set(ROCM_DEFAULT OFF)
endif()
if(KOKKOS_ENABLE_HIP)
  set(ROCTHRUST_DEFAULT ON)
else()
  set(ROCTHRUST_DEFAULT OFF)
endif()
kokkos_tpl_option(ROCM ${ROCM_DEFAULT})
kokkos_tpl_option(ROCTHRUST ${ROCTHRUST_DEFAULT})
if(Kokkos_ENABLE_ROCTHRUST)
  include(CheckCXXSourceCompiles)
  check_cxx_source_compiles(
    "
    #include <ios>
    int main() {
      static_assert(_GLIBCXX_RELEASE < 9);
      return 0;
    }
    "
    Kokkos_ENABLE_IMPL_SKIP_NO_RTTI_FLAG
  )
endif()

if(KOKKOS_ENABLE_SYCL)
  set(ONEDPL_DEFAULT ON)
else()
  set(ONEDPL_DEFAULT OFF)
endif()
kokkos_tpl_option(ONEDPL ${ONEDPL_DEFAULT})

if(WIN32)
  set(LIBDL_DEFAULT Off)
else()
  set(LIBDL_DEFAULT On)
endif()
kokkos_tpl_option(LIBDL ${LIBDL_DEFAULT} TRIBITS DLlib)

if(Trilinos_ENABLE_Kokkos AND TPL_ENABLE_HPX)
  set(HPX_DEFAULT ON)
else()
  set(HPX_DEFAULT OFF)
endif()
kokkos_tpl_option(HPX ${HPX_DEFAULT})

kokkos_tpl_option(THREADS ${Kokkos_ENABLE_THREADS} TRIBITS Pthread)

if(Trilinos_ENABLE_Kokkos AND TPL_ENABLE_quadmath)
  set(LIBQUADMATH_DEFAULT ON)
else()
  set(LIBQUADMATH_DEFAULT OFF)
endif()
kokkos_tpl_option(LIBQUADMATH ${LIBQUADMATH_DEFAULT} TRIBITS quadmath)

#Make sure we use our local FindKokkosCuda.cmake
kokkos_import_tpl(HPX INTERFACE)
kokkos_import_tpl(CUDA INTERFACE)
kokkos_import_tpl(HWLOC)
kokkos_import_tpl(LIBDL)
if(NOT WIN32)
  kokkos_import_tpl(THREADS INTERFACE)
endif()
if(NOT KOKKOS_ENABLE_COMPILE_AS_CMAKE_LANGUAGE)
  kokkos_import_tpl(ROCM INTERFACE)
endif()
kokkos_import_tpl(ONEDPL INTERFACE)
kokkos_import_tpl(LIBQUADMATH)
kokkos_import_tpl(ROCTHRUST)

if(Kokkos_ENABLE_DESUL_ATOMICS_EXTERNAL)
  find_package(desul REQUIRED COMPONENTS atomics)
  kokkos_export_cmake_tpl(desul REQUIRED COMPONENTS atomics)
endif()

if(Kokkos_ENABLE_IMPL_MDSPAN AND Kokkos_ENABLE_MDSPAN_EXTERNAL)
  find_package(mdspan REQUIRED)
  kokkos_export_cmake_tpl(mdspan REQUIRED)
endif()

if(Kokkos_ENABLE_OPENMP)
  find_package(OpenMP 3.0 REQUIRED COMPONENTS CXX)
  kokkos_export_cmake_tpl(OpenMP REQUIRED COMPONENTS CXX)
  if(Kokkos_ENABLE_HIP AND KOKKOS_COMPILE_LANGUAGE STREQUAL HIP)
    global_append(KOKKOS_AMDGPU_OPTIONS ${OpenMP_CXX_FLAGS})
  endif()
  if(Kokkos_ENABLE_CUDA AND KOKKOS_COMPILE_LANGUAGE STREQUAL CUDA)
    global_append(KOKKOS_CUDA_OPTIONS -Xcompiler ${OpenMP_CXX_FLAGS})
  endif()
endif()

#Convert list to newlines (which CMake doesn't always like in cache variables)
string(REPLACE ";" "\n" KOKKOS_TPL_EXPORT_TEMP "${KOKKOS_TPL_EXPORTS}")
#Convert to a regular variable
unset(KOKKOS_TPL_EXPORTS CACHE)
set(KOKKOS_TPL_EXPORTS ${KOKKOS_TPL_EXPORT_TEMP})
