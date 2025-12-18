#I have to leave these here for tribits
kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkos_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

file(GLOB SIMD_HEADERS *.hpp)
file(GLOB SIMD_SOURCES *.cpp)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION ${KOKKOS_HEADER_DIR}
  FILES_MATCHING
  PATTERN "*.hpp"
)

set(SIMD_MODULES Kokkos_SIMD.cppm Kokkos_SIMD_Impl.cppm)

#-----------------------------------------------------------------------------

# We have to pass the sources in here for Tribits
# These will get ignored for standalone CMake and a true interface library made
kokkos_add_library(
  kokkossimd
  STATIC # MSVC can't deal with empty shared libraries
  SOURCES
  ${SIMD_SOURCES}
  HEADERS
  ${SIMD_HEADERS}
  MODULE_INTERFACE
  ${SIMD_MODULES}
)
kokkos_lib_include_directories(
  kokkossimd ${KOKKOS_TOP_BUILD_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
)

kokkos_link_internal_library(kokkossimd kokkoscore)

# If ARM-SVE, check and warn for available arm_neon_sve_bridge.h
if(KOKKOS_ARCH_ARM_SVE)
  include(CheckIncludeFileCXX)
  check_include_file_cxx(arm_neon_sve_bridge.h KOKKOS_COMPILER_SUPPORTS_ARM_NEON_SVE_BRIDGE)
  if(NOT KOKKOS_COMPILER_SUPPORTS_ARM_NEON_SVE_BRIDGE)
    message(
      WARNING "You are using a compiler without NEON-SVE bridge header "
              "(arm_neon_sve_bridge.h). Kokkos will use its own "
              "implementation of these functions, which could be slower "
              "than native ones from a more recent SVE compiler. It is "
              "recommended to upgrade your compiler.\n"
              "NOTE: It was also observed that, in some LLVM versions "
              "(14, 15), despite the macro __ARM_NEON_SVE_BRIDGE was "
              "defined, the header arm_neon_sve_bridge.h is still not "
              "available (see https://godbolt.org/z/7jqWvzvWY )."
    )
  endif()
endif()
