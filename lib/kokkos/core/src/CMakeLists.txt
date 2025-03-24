kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${KOKKOS_TOP_BUILD_DIR})
if(NOT desul_FOUND)
  if(KOKKOS_ENABLE_CUDA)
    set(DESUL_ATOMICS_ENABLE_CUDA ON)
  endif()
  if(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
    set(DESUL_ATOMICS_ENABLE_CUDA_SEPARABLE_COMPILATION ON)
  endif()
  if(KOKKOS_ENABLE_HIP)
    set(DESUL_ATOMICS_ENABLE_HIP ON)
  endif()
  if(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
    set(DESUL_ATOMICS_ENABLE_HIP_SEPARABLE_COMPILATION ON)
  endif()
  if(KOKKOS_ENABLE_SYCL)
    set(DESUL_ATOMICS_ENABLE_SYCL ON)
    if(KOKKOS_IMPL_SYCL_DEVICE_GLOBAL_SUPPORTED AND NOT KOKKOS_IMPL_HAVE_SYCL_EXT_ONEAPI_DEVICE_GLOBAL)
      set(DESUL_ATOMICS_ENABLE_SYCL_SEPARABLE_COMPILATION ON)
    endif()
  endif()
  if(KOKKOS_ENABLE_OPENMPTARGET)
    set(DESUL_ATOMICS_ENABLE_OPENMP ON) # not a typo Kokkos OpenMPTarget -> Desul OpenMP
  endif()
  if(KOKKOS_ENABLE_OPENACC)
    # FIXME_OPENACC FIXME_CLACC - Below condition will be removed if Clacc can compile atomics.
    if(KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
      set(DESUL_ATOMICS_ENABLE_OPENACC ON)
    endif()
  endif()
  configure_file(
    ${KOKKOS_SOURCE_DIR}/tpls/desul/Config.hpp.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/desul/atomics/Config.hpp
  )
  kokkos_include_directories(${KOKKOS_SOURCE_DIR}/tpls/desul/include)
endif()

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION ${KOKKOS_HEADER_DIR}
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h"
)

set(KOKKOS_CORE_SRCS)
append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/impl/*.cpp)
set(KOKKOS_CORE_HEADERS)
append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/impl/*.hpp)

if(KOKKOS_ENABLE_CUDA)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Cuda/*.cpp)
  if(NOT Kokkos_ENABLE_DEPRECATED_CODE_4)
    list(REMOVE_ITEM KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Cuda/Kokkos_Cuda_Task.cpp)
  endif()
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Cuda/*.hpp)
endif()

if(KOKKOS_ENABLE_OPENMP)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMP/*.cpp)
  if(NOT Kokkos_ENABLE_DEPRECATED_CODE_4)
    list(REMOVE_ITEM KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMP/Kokkos_OpenMP_Task.cpp)
  endif()
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMP/*.hpp)
endif()

if(KOKKOS_ENABLE_OPENMPTARGET)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMPTarget/*.cpp)
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMPTarget/*.hpp)
endif()

if(KOKKOS_ENABLE_OPENACC)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/OpenACC/*.cpp)
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/OpenACC/*.hpp)
endif()

if(KOKKOS_ENABLE_THREADS)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Threads/*.cpp)
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Threads/*.hpp)
endif()

if(KOKKOS_ENABLE_HIP)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/HIP/*.cpp)
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/HIP/*.hpp)
endif()

if(KOKKOS_ENABLE_HPX)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/HPX/*.cpp)
  if(NOT Kokkos_ENABLE_DEPRECATED_CODE_4)
    list(REMOVE_ITEM KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/HPX/Kokkos_HPX_Task.cpp)
  endif()
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/HPX/*.hpp)
endif()

if(KOKKOS_ENABLE_SERIAL)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Serial/*.cpp)
  if(NOT Kokkos_ENABLE_DEPRECATED_CODE_4)
    list(REMOVE_ITEM KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Serial/Kokkos_Serial_Task.cpp)
  endif()
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Serial/*.hpp)
endif()

if(KOKKOS_ENABLE_SYCL)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/SYCL/*.cpp)
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/SYCL/*.hpp)
endif()

if(NOT desul_FOUND)
  if(KOKKOS_ENABLE_CUDA)
    append_glob(KOKKOS_CORE_SRCS ${KOKKOS_SOURCE_DIR}/tpls/desul/src/Lock_Array_CUDA.cpp)
  elseif(KOKKOS_ENABLE_HIP)
    append_glob(KOKKOS_CORE_SRCS ${KOKKOS_SOURCE_DIR}/tpls/desul/src/Lock_Array_HIP.cpp)
  elseif(KOKKOS_ENABLE_SYCL)
    append_glob(KOKKOS_CORE_SRCS ${KOKKOS_SOURCE_DIR}/tpls/desul/src/Lock_Array_SYCL.cpp)
  endif()
  append_glob(KOKKOS_CORE_HEADERS ${KOKKOS_SOURCE_DIR}/tpls/desul/include/desul/*.hpp)
  append_glob(KOKKOS_CORE_HEADERS ${KOKKOS_SOURCE_DIR}/tpls/desul/include/desul/*/*.hpp)
  append_glob(KOKKOS_CORE_HEADERS ${KOKKOS_SOURCE_DIR}/tpls/desul/include/desul/*/*/*.hpp)
  append_glob(KOKKOS_CORE_HEADERS ${KOKKOS_SOURCE_DIR}/tpls/desul/include/*/*/*.inc*)
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/desul/*.hpp)

  install(
    DIRECTORY "${KOKKOS_SOURCE_DIR}/tpls/desul/include/desul" "${CMAKE_CURRENT_BINARY_DIR}/desul"
    DESTINATION ${KOKKOS_HEADER_DIR}
    FILES_MATCHING
    PATTERN "*.inc"
    PATTERN "*.inc_*"
    PATTERN "*.hpp"
  )

  message(STATUS "Using internal desul_atomics copy")
else()
  message(STATUS "Using external desul_atomics install found at:")
  message(STATUS "  " ${desul_DIR})
endif()

kokkos_add_library(
  kokkoscore SOURCES ${KOKKOS_CORE_SRCS} HEADERS ${KOKKOS_CORE_HEADERS}
  ADD_BUILD_OPTIONS # core should be given all the necessary compiler/linker flags
)

kokkos_lib_include_directories(
  kokkoscore ${KOKKOS_TOP_BUILD_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
)
if(NOT desul_FOUND)
  kokkos_lib_include_directories(kokkoscore ${KOKKOS_SOURCE_DIR}/tpls/desul/include)
endif()

if(Kokkos_ENABLE_IMPL_MDSPAN)
  message(STATUS "Experimental mdspan support is enabled")

  # Some compilers now include mdspan... we just flag on their version
  # for now until we can get some compiler detection support
  include(CheckIncludeFileCXX)
  check_include_file_cxx(experimental/mdspan KOKKOS_COMPILER_SUPPORTS_EXPERIMENTAL_MDSPAN)
  check_include_file_cxx(mdspan KOKKOS_COMPILER_SUPPORTS_MDSPAN)

  if(Kokkos_ENABLE_MDSPAN_EXTERNAL)
    message(STATUS "Using external mdspan")
    target_link_libraries(kokkoscore PUBLIC std::mdspan)
  elseif(KOKKOS_COMPILER_SUPPORTS_MDSPAN AND NOT Kokkos_ENABLE_IMPL_SKIP_COMPILER_MDSPAN)
    message(STATUS "Using compiler-supplied mdspan")
  elseif(KOKKOS_COMPILER_SUPPORTS_EXPERIMENTAL_MDSPAN AND NOT Kokkos_ENABLE_IMPL_SKIP_COMPILER_MDSPAN)
    message(STATUS "Using compiler-supplied experimental/mdspan")
  else()
    kokkos_lib_include_directories(kokkoscore ${KOKKOS_SOURCE_DIR}/tpls/mdspan/include)

    append_glob(KOKKOS_CORE_HEADERS ${KOKKOS_SOURCE_DIR}/tpls/mdspan/include/experimental/__p0009_bits/*.hpp)
    append_glob(KOKKOS_CORE_HEADERS ${KOKKOS_SOURCE_DIR}/tpls/mdspan/include/experimental/mdspan)

    install(
      DIRECTORY "${KOKKOS_SOURCE_DIR}/tpls/mdspan/include/"
      DESTINATION ${KOKKOS_HEADER_DIR}
      FILES_MATCHING
      PATTERN "mdspan"
      PATTERN "*.hpp"
    )
    message(STATUS "Using internal mdspan directory ${KOKKOS_SOURCE_DIR}/tpls/mdspan/include")
  endif()
endif()

kokkos_link_tpl(kokkoscore PUBLIC HWLOC)
kokkos_link_tpl(kokkoscore PUBLIC CUDA)
kokkos_link_tpl(kokkoscore PUBLIC HPX)
kokkos_link_tpl(kokkoscore PUBLIC LIBDL)
# On *nix-like systems (Linux, macOS) we need pthread for C++ std::thread
if(NOT WIN32)
  kokkos_link_tpl(kokkoscore PUBLIC THREADS)
endif()
if(NOT KOKKOS_ENABLE_COMPILE_AS_CMAKE_LANGUAGE)
  kokkos_link_tpl(kokkoscore PUBLIC ROCM)
endif()

# FIXME: We need a proper solution to figure out whether to enable
#        libatomic
# Most compilers only require libatomic for 128-bit CAS
# I (CT) had removed 128bit CAS from desul to not need libatomic.
if(KOKKOS_ENABLE_OPENMPTARGET)
  target_link_libraries(kokkoscore PUBLIC atomic)
endif()

if(desul_FOUND)
  target_link_libraries(kokkoscore PUBLIC desul_atomics)
endif()

if(Kokkos_ENABLE_OPENMP)
  target_link_libraries(kokkoscore PUBLIC OpenMP::OpenMP_CXX)
endif()

kokkos_link_tpl(kokkoscore PUBLIC LIBQUADMATH)
