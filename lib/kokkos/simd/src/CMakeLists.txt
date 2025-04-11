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

#-----------------------------------------------------------------------------

# We have to pass the sources in here for Tribits
# These will get ignored for standalone CMake and a true interface library made
kokkos_add_library(kokkossimd SOURCES ${SIMD_SOURCES} HEADERS ${SIMD_HEADERS})
kokkos_lib_include_directories(
  kokkossimd ${KOKKOS_TOP_BUILD_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
)
