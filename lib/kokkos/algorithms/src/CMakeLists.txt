#I have to leave these here for tribits
kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkos_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

file(GLOB ALGO_HEADERS *.hpp)
file(GLOB ALGO_SOURCES *.cpp)
append_glob(ALGO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/std_algorithms/*.hpp)
append_glob(ALGO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/std_algorithms/impl/*.hpp)

set(ALGO_MODULES Kokkos_Random.cppm Kokkos_Sort.cppm Kokkos_StdAlgorithms.cppm)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION ${KOKKOS_HEADER_DIR}
  FILES_MATCHING
  PATTERN "*.hpp"
)

#-----------------------------------------------------------------------------

kokkos_add_library(
  kokkosalgorithms
  HEADERS
  ${ALGO_HEADERS}
  MODULE_INTERFACE
  ${ALGO_MODULES}
  SOURCES
  ${ALGO_SOURCES}
)
kokkos_lib_include_directories(
  kokkosalgorithms ${KOKKOS_TOP_BUILD_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(kokkosalgorithms PUBLIC kokkoscore)

kokkos_link_tpl(kokkosalgorithms PUBLIC ROCTHRUST)
kokkos_link_tpl(kokkosalgorithms PUBLIC ONEDPL)
