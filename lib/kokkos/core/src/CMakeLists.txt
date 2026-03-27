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
  if(BUILD_SHARED_LIBS)
    set(DESUL_IMPL_BUILD_SHARED_LIBS ON)
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
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Cuda/*.hpp)
endif()

if(KOKKOS_ENABLE_OPENMP)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMP/*.cpp)
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
  append_glob(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/HPX/*.hpp)
endif()

if(KOKKOS_ENABLE_SERIAL)
  append_glob(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Serial/*.cpp)
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

  file(STRINGS "${KOKKOS_SOURCE_DIR}/tpls/desul-hash.txt" KOKKOS_DESUL_HASH)
  global_set(KOKKOS_DESUL_VERSION "${KOKKOS_DESUL_HASH}")
  message(STATUS "Using bundled desul_atomics copy (desul/desul@${KOKKOS_DESUL_HASH})")
else()
  global_set(KOKKOS_DESUL_VERSION "unknown")
  message(STATUS "Using external desul_atomics install found at:")
  message(STATUS "  " ${desul_DIR})
endif()

if(Kokkos_ENABLE_EXPERIMENTAL_CXX20_MODULES)
  set(KOKKOS_CORE_MODULE_FILES Kokkos_Core.cppm Kokkos_Core_Impl.cppm)
endif()

kokkos_add_library(
  kokkoscore
  SOURCES
  ${KOKKOS_CORE_SRCS}
  MODULE_INTERFACE
  ${KOKKOS_CORE_MODULE_FILES}
  HEADERS
  ${KOKKOS_CORE_HEADERS}
  ADD_BUILD_OPTIONS # core should be given all the necessary compiler/linker flags
)

kokkos_lib_include_directories(
  kokkoscore ${KOKKOS_TOP_BUILD_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
)
if(NOT desul_FOUND)
  target_include_directories(kokkoscore SYSTEM PUBLIC $<BUILD_INTERFACE:${KOKKOS_SOURCE_DIR}/tpls/desul/include>)
  if(BUILD_SHARED_LIBS)
    target_compile_definitions(kokkoscore PRIVATE DESUL_IMPL_EXPORT_SYMBOLS)
  endif()
endif()

if(Kokkos_ENABLE_IMPL_MDSPAN)
  if(Kokkos_ENABLE_MDSPAN_EXTERNAL)
    global_set(KOKKOS_MDSPAN_VERSION "unknown")
    message(STATUS "Using external mdspan")
    target_link_libraries(kokkoscore PUBLIC std::mdspan)
  else()
    target_include_directories(kokkoscore SYSTEM PUBLIC $<BUILD_INTERFACE:${KOKKOS_SOURCE_DIR}/tpls/mdspan/include>)

    append_glob(KOKKOS_CORE_HEADERS ${KOKKOS_SOURCE_DIR}/tpls/mdspan/include/experimental/__p0009_bits/*.hpp)
    append_glob(KOKKOS_CORE_HEADERS ${KOKKOS_SOURCE_DIR}/tpls/mdspan/include/experimental/mdspan)

    install(
      DIRECTORY "${KOKKOS_SOURCE_DIR}/tpls/mdspan/include/"
      DESTINATION ${KOKKOS_HEADER_DIR}
      FILES_MATCHING
      PATTERN "mdspan"
      PATTERN "*.hpp"
    )

    file(STRINGS "${KOKKOS_SOURCE_DIR}/tpls/mdspan-hash.txt" KOKKOS_MDSPAN_HASH)
    global_set(KOKKOS_MDSPAN_VERSION "${KOKKOS_MDSPAN_HASH}")
    message(STATUS "Using bundled mdspan copy (kokkos/mdspan@${KOKKOS_MDSPAN_HASH})")
  endif()
endif()

if(Kokkos_ENABLE_IMPL_VIEW_LEGACY)
  global_set(KOKKOS_MDSPAN_VERSION "not applicable")
endif()

kokkos_link_tpl(kokkoscore PUBLIC HWLOC)
kokkos_link_tpl(kokkoscore PUBLIC CUDA)
kokkos_link_tpl(kokkoscore PUBLIC HPX)
# On *nix-like systems (Linux, macOS) we need pthread for C++ std::thread
if(NOT WIN32)
  kokkos_link_tpl(kokkoscore PUBLIC THREADS)
endif()
if(NOT KOKKOS_ENABLE_COMPILE_AS_CMAKE_LANGUAGE)
  kokkos_link_tpl(kokkoscore PUBLIC ROCM)
endif()
kokkos_link_tpl(kokkoscore PUBLIC LIBQUADMATH)

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

if(Kokkos_ENABLE_LIBDL)
  target_link_libraries(kokkoscore PUBLIC ${CMAKE_DL_LIBS})
endif()
