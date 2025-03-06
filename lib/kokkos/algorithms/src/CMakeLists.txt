#I have to leave these here for tribits
kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkos_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

file(GLOB ALGO_HEADERS *.hpp)
file(GLOB ALGO_SOURCES *.cpp)
append_glob(ALGO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/std_algorithms/*.hpp)
append_glob(ALGO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/std_algorithms/impl/*.hpp)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION ${KOKKOS_HEADER_DIR}
  FILES_MATCHING
  PATTERN "*.hpp"
)

#-----------------------------------------------------------------------------

# We have to pass the sources in here for Tribits
# These will get ignored for standalone CMake and a true interface library made
kokkos_add_interface_library(kokkosalgorithms NOINSTALLHEADERS ${ALGO_HEADERS} SOURCES ${ALGO_SOURCES})
kokkos_lib_include_directories(
  kokkosalgorithms ${KOKKOS_TOP_BUILD_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
)

kokkos_link_tpl(kokkoscontainers PUBLIC ROCTHRUST)
kokkos_link_tpl(kokkoscore PUBLIC ONEDPL)
