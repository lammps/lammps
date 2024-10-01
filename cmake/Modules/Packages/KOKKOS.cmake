########################################################################
# As of version 4.0.0 Kokkos requires C++17
if(CMAKE_CXX_STANDARD LESS 17)
  message(FATAL_ERROR "The KOKKOS package requires the C++ standard to
be set to at least C++17")
endif()

########################################################################
# consistency checks and Kokkos options/settings required by LAMMPS
if(Kokkos_ENABLE_CUDA)
  option(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC "CUDA asynchronous malloc support" OFF)
  mark_as_advanced(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC)
  if(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC)
    message(STATUS "KOKKOS: CUDA malloc async support enabled")
  else()
    message(STATUS "KOKKOS: CUDA malloc async support disabled")
  endif()
endif()
if(Kokkos_ENABLE_HIP)
  option(Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS "Enable multiple kernel instantiations with HIP" ON)
  mark_as_advanced(Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS)
  option(Kokkos_ENABLE_ROCTHRUST "Use RoCThrust library" ON)
  mark_as_advanced(Kokkos_ENABLE_ROCTHRUST)

  if(Kokkos_ARCH_AMD_GFX942 OR Kokkos_ARCH_AMD_GFX940)
    option(Kokkos_ENABLE_IMPL_HIP_UNIFIED_MEMORY "Enable unified memory with HIP" ON)
    mark_as_advanced(Kokkos_ENABLE_IMPL_HIP_UNIFIED_MEMORY)
  endif()
endif()
# Adding OpenMP compiler flags without the checks done for
# BUILD_OMP can result in compile failures. Enforce consistency.
if(Kokkos_ENABLE_OPENMP)
  if(NOT BUILD_OMP)
    message(FATAL_ERROR "Must enable BUILD_OMP with Kokkos_ENABLE_OPENMP")
  endif()
endif()

if(Kokkos_ENABLE_SERIAL)
  if(NOT (Kokkos_ENABLE_OPENMP OR Kokkos_ENABLE_THREADS OR
    Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL
    OR Kokkos_ENABLE_OPENMPTARGET))
  option(Kokkos_ENABLE_ATOMICS_BYPASS "Disable atomics for Kokkos Serial Backend" ON)
  mark_as_advanced(Kokkos_ENABLE_ATOMICS_BYPASS)
  endif()
endif()
########################################################################

option(EXTERNAL_KOKKOS "Build against external kokkos library" OFF)
option(DOWNLOAD_KOKKOS "Download the KOKKOS library instead of using the bundled one" OFF)
if(DOWNLOAD_KOKKOS)
  # extract Kokkos-related variables and values so we can forward them to the Kokkos library build
  get_cmake_property(_VARS VARIABLES)
  list(FILTER _VARS INCLUDE REGEX ^Kokkos_)
  foreach(_VAR IN LISTS _VARS)
    list(APPEND KOKKOS_LIB_BUILD_ARGS "-D${_VAR}=${${_VAR}}")
  endforeach()
  message(STATUS "KOKKOS download requested - we will build our own")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
  # build KOKKOS downloaded libraries as static libraries but with PIC, if needed
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DBUILD_SHARED_LIBS=OFF")
  if(CMAKE_REQUEST_PIC)
    list(APPEND KOKKOS_LIB_BUILD_ARGS ${CMAKE_REQUEST_PIC})
  endif()
  # append other CMake variables that need to be forwarded to CMAKE_ARGS
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_INSTALL_LIBDIR=lib")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
  include(ExternalProject)
  set(KOKKOS_URL "https://github.com/kokkos/kokkos/archive/4.4.01.tar.gz" CACHE STRING "URL for KOKKOS tarball")
  set(KOKKOS_MD5 "de6ee80d00b6212b02bfb7f1e71a8392" CACHE STRING "MD5 checksum of KOKKOS tarball")
  mark_as_advanced(KOKKOS_URL)
  mark_as_advanced(KOKKOS_MD5)
  GetFallbackURL(KOKKOS_URL KOKKOS_FALLBACK)

  ExternalProject_Add(kokkos_build
    URL     ${KOKKOS_URL} ${KOKKOS_FALLBACK}
    URL_MD5 ${KOKKOS_MD5}
    CMAKE_ARGS ${KOKKOS_LIB_BUILD_ARGS}
    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/libkokkoscore.a <INSTALL_DIR>/lib/libkokkoscontainers.a
  )
  ExternalProject_get_property(kokkos_build INSTALL_DIR)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
  add_library(LAMMPS::KOKKOSCORE UNKNOWN IMPORTED)
  add_library(LAMMPS::KOKKOSCONTAINERS UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::KOKKOSCORE PROPERTIES
    IMPORTED_LOCATION "${INSTALL_DIR}/lib/libkokkoscore.a"
    INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/include"
    INTERFACE_LINK_LIBRARIES ${CMAKE_DL_LIBS})
  set_target_properties(LAMMPS::KOKKOSCONTAINERS PROPERTIES
    IMPORTED_LOCATION "${INSTALL_DIR}/lib/libkokkoscontainers.a")
  target_link_libraries(lammps PRIVATE LAMMPS::KOKKOSCORE LAMMPS::KOKKOSCONTAINERS)
  add_dependencies(LAMMPS::KOKKOSCORE kokkos_build)
  add_dependencies(LAMMPS::KOKKOSCONTAINERS kokkos_build)
elseif(EXTERNAL_KOKKOS)
  find_package(Kokkos 4.4.01 REQUIRED CONFIG)
  target_link_libraries(lammps PRIVATE Kokkos::kokkos)
else()
  set(LAMMPS_LIB_KOKKOS_SRC_DIR ${LAMMPS_LIB_SOURCE_DIR}/kokkos)
  set(LAMMPS_LIB_KOKKOS_BIN_DIR ${LAMMPS_LIB_BINARY_DIR}/kokkos)
  # build KOKKOS internal libraries as static libraries but with PIC, if needed
  if(BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS_WAS_ON YES)
    set(BUILD_SHARED_LIBS OFF)
  endif()
  if(CMAKE_REQUEST_PIC)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  endif()
  add_subdirectory(${LAMMPS_LIB_KOKKOS_SRC_DIR} ${LAMMPS_LIB_KOKKOS_BIN_DIR} EXCLUDE_FROM_ALL)

  set(Kokkos_INCLUDE_DIRS ${LAMMPS_LIB_KOKKOS_SRC_DIR}/core/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/containers/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/algorithms/src
                          ${LAMMPS_LIB_KOKKOS_BIN_DIR})
  target_include_directories(lammps PRIVATE ${Kokkos_INCLUDE_DIRS})
  target_link_libraries(lammps PRIVATE kokkos)
  if(BUILD_SHARED_LIBS_WAS_ON)
    set(BUILD_SHARED_LIBS ON)
  endif()
endif()
target_compile_definitions(lammps PUBLIC $<BUILD_INTERFACE:LMP_KOKKOS>)

set(KOKKOS_PKG_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/KOKKOS)
set(KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/atom_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/atom_map_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/atom_vec_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/comm_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/comm_tiled_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/min_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/min_linesearch_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/neighbor_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/neigh_list_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/neigh_bond_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/fix_nh_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/nbin_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/npair_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/npair_halffull_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/domain_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/modify_kokkos.cpp)

# fix wall/gran has been refactored in an incompatible way. Use old version of base class for now
if(PKG_GRANULAR)
  list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/fix_wall_gran_old.cpp)
endif()

if(PKG_KSPACE)
  list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/fft3d_kokkos.cpp
                                 ${KOKKOS_PKG_SOURCES_DIR}/grid3d_kokkos.cpp
                                 ${KOKKOS_PKG_SOURCES_DIR}/remap_kokkos.cpp)
  set(FFT_KOKKOS "KISS" CACHE STRING "FFT library for Kokkos-enabled KSPACE package")
  set(FFT_KOKKOS_VALUES KISS FFTW3 MKL NVPL HIPFFT CUFFT MKL_GPU)
  set_property(CACHE FFT_KOKKOS PROPERTY STRINGS ${FFT_KOKKOS_VALUES})
  validate_option(FFT_KOKKOS FFT_KOKKOS_VALUES)
  string(TOUPPER ${FFT_KOKKOS} FFT_KOKKOS)

  if(Kokkos_ENABLE_CUDA)
    if(NOT ((FFT_KOKKOS STREQUAL "KISS") OR (FFT_KOKKOS STREQUAL "CUFFT")))
      message(FATAL_ERROR "The CUDA backend of Kokkos requires either KISS FFT or CUFFT.")
    elseif(FFT_KOKKOS STREQUAL "KISS")
      message(WARNING "Using KISS FFT with the CUDA backend of Kokkos may be sub-optimal.")
    elseif(FFT_KOKKOS STREQUAL "CUFFT")
      find_package(CUDAToolkit REQUIRED)
      target_link_libraries(lammps PRIVATE CUDA::cufft)
    endif()
  elseif(Kokkos_ENABLE_HIP)
    if(NOT ((FFT_KOKKOS STREQUAL "KISS") OR (FFT_KOKKOS STREQUAL "HIPFFT")))
      message(FATAL_ERROR "The HIP backend of Kokkos requires either KISS FFT or HIPFFT.")
    elseif(FFT_KOKKOS STREQUAL "KISS")
      message(WARNING "Using KISS FFT with the HIP backend of Kokkos may be sub-optimal.")
      target_compile_definitions(lammps PRIVATE -DFFT_KOKKOS_KISS)
    elseif(FFT_KOKKOS STREQUAL "HIPFFT")
      include(DetectHIPInstallation)
      find_package(hipfft REQUIRED)
      target_link_libraries(lammps PRIVATE hip::hipfft)
    endif()
  elseif(FFT_KOKKOS STREQUAL "MKL_GPU")
    if(NOT Kokkos_ENABLE_SYCL)
      message(FATAL_ERROR "Using MKL_GPU FFT currently requires the SYCL backend of Kokkos.")
    endif()
    find_package(MKL REQUIRED)
    target_link_libraries(lammps PRIVATE mkl_sycl_dft mkl_intel_ilp64 mkl_tbb_thread mkl_core tbb)
  elseif(FFT_KOKKOS STREQUAL "MKL")
    find_package(MKL REQUIRED)
  elseif(FFT_KOKKOS STREQUAL "NVPL")
      find_package(nvpl_fft REQUIRED)
      target_link_libraries(lammps PRIVATE nvpl::fftw)
  endif()
  target_compile_definitions(lammps PRIVATE -DFFT_KOKKOS_${FFT_KOKKOS})
endif()

if(PKG_ML-IAP)
  list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/mliap_data_kokkos.cpp
                                 ${KOKKOS_PKG_SOURCES_DIR}/mliap_descriptor_so3_kokkos.cpp
                                 ${KOKKOS_PKG_SOURCES_DIR}/mliap_model_linear_kokkos.cpp
                                 ${KOKKOS_PKG_SOURCES_DIR}/mliap_model_python_kokkos.cpp
                                 ${KOKKOS_PKG_SOURCES_DIR}/mliap_unified_kokkos.cpp
                                 ${KOKKOS_PKG_SOURCES_DIR}/mliap_so3_kokkos.cpp)

  # Add KOKKOS version of ML-IAP Python coupling if non-KOKKOS version is included
  if(MLIAP_ENABLE_PYTHON AND Cythonize_EXECUTABLE)
    file(GLOB MLIAP_KOKKOS_CYTHON_SRC CONFIGURE_DEPENDS ${LAMMPS_SOURCE_DIR}/KOKKOS/*.pyx)
    foreach(MLIAP_CYTHON_FILE ${MLIAP_KOKKOS_CYTHON_SRC})
      get_filename_component(MLIAP_CYTHON_BASE ${MLIAP_CYTHON_FILE} NAME_WE)
      add_custom_command(OUTPUT  ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.cpp ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.h
              COMMAND            ${CMAKE_COMMAND} -E copy_if_different ${MLIAP_CYTHON_FILE} ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.pyx
              COMMAND            ${Cythonize_EXECUTABLE} -3 ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.pyx
              WORKING_DIRECTORY  ${MLIAP_BINARY_DIR}
              MAIN_DEPENDENCY    ${MLIAP_CYTHON_FILE}
              COMMENT "Generating C++ sources with cythonize...")
      list(APPEND KOKKOS_PKG_SOURCES ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.cpp)
    endforeach()
  endif()
endif()

if(PKG_PHONON)
  list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/dynamical_matrix_kokkos.cpp)
  list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/third_order_kokkos.cpp)
endif()

set_property(GLOBAL PROPERTY "KOKKOS_PKG_SOURCES" "${KOKKOS_PKG_SOURCES}")

# detects styles which have KOKKOS version
RegisterStylesExt(${KOKKOS_PKG_SOURCES_DIR} kokkos KOKKOS_PKG_SOURCES)

# register kokkos-only styles
RegisterNBinStyle(${KOKKOS_PKG_SOURCES_DIR}/nbin_kokkos.h)
RegisterNPairStyle(${KOKKOS_PKG_SOURCES_DIR}/npair_kokkos.h)
RegisterNPairStyle(${KOKKOS_PKG_SOURCES_DIR}/npair_halffull_kokkos.h)

if(PKG_DPD-REACT)
  get_property(KOKKOS_PKG_SOURCES GLOBAL PROPERTY KOKKOS_PKG_SOURCES)
  list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/npair_ssa_kokkos.cpp)
  RegisterNPairStyle(${KOKKOS_PKG_SOURCES_DIR}/npair_ssa_kokkos.h)
  set_property(GLOBAL PROPERTY "KOKKOS_PKG_SOURCES" "${KOKKOS_PKG_SOURCES}")
endif()

get_property(KOKKOS_PKG_SOURCES GLOBAL PROPERTY KOKKOS_PKG_SOURCES)

target_sources(lammps PRIVATE ${KOKKOS_PKG_SOURCES})
target_include_directories(lammps PUBLIC $<BUILD_INTERFACE:${KOKKOS_PKG_SOURCES_DIR}>)
