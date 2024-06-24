#I have to leave these here for tribits
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

FILE(GLOB ALGO_HEADERS *.hpp)
FILE(GLOB ALGO_SOURCES *.cpp)
APPEND_GLOB(ALGO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/std_algorithms/*.hpp)
APPEND_GLOB(ALGO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/std_algorithms/impl/*.hpp)

INSTALL (
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION ${KOKKOS_HEADER_DIR}
  FILES_MATCHING PATTERN "*.hpp"
)

#-----------------------------------------------------------------------------

# We have to pass the sources in here for Tribits
# These will get ignored for standalone CMake and a true interface library made
KOKKOS_ADD_INTERFACE_LIBRARY(
  kokkosalgorithms
  NOINSTALLHEADERS ${ALGO_HEADERS}
  SOURCES ${ALGO_SOURCES}
)
KOKKOS_LIB_INCLUDE_DIRECTORIES(kokkosalgorithms
  ${KOKKOS_TOP_BUILD_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

KOKKOS_LINK_TPL(kokkoscontainers PUBLIC ROCTHRUST)
KOKKOS_LINK_TPL(kokkoscore PUBLIC ONEDPL)
