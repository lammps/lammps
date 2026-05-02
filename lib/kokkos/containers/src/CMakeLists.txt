#need these here for now
kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkos_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

set(KOKKOS_CONTAINERS_SRCS)
append_glob(KOKKOS_CONTAINERS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/impl/*.cpp)
set(KOKKOS_CONTAINER_HEADERS)
append_glob(KOKKOS_CONTAINERS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/impl/*.hpp)
append_glob(KOKKOS_CONTAINERS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)

set(KOKKOS_CONTAINERS_MODULES
    Kokkos_Bitset.cppm
    Kokkos_DualView.cppm
    Kokkos_DynamicView.cppm
    Kokkos_DynRankView.cppm
    Kokkos_DynRankView_Impl.cppm
    Kokkos_ErrorReporter.cppm
    Kokkos_Functional.cppm
    Kokkos_OffsetView.cppm
    Kokkos_ScatterView.cppm
    Kokkos_ScatterView_Impl.cppm
    Kokkos_UnorderedMap.cppm
    Kokkos_UnorderedMap_Impl.cppm
)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION ${KOKKOS_HEADER_DIR}
  FILES_MATCHING
  PATTERN "*.hpp"
)

kokkos_add_library(
  kokkoscontainers
  SOURCES
  ${KOKKOS_CONTAINERS_SRCS}
  HEADERS
  ${KOKKOS_CONTAINERS_HEADERS}
  MODULE_INTERFACE
  ${KOKKOS_CONTAINERS_MODULES}
)

kokkos_lib_include_directories(
  kokkoscontainers ${KOKKOS_TOP_BUILD_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
)
kokkos_link_internal_library(kokkoscontainers kokkoscore)

#-----------------------------------------------------------------------------
