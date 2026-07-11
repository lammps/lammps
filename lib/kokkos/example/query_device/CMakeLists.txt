kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkos_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES "")
file(GLOB SOURCES *.cpp)

kokkos_add_executable(query_device SOURCES ${SOURCES})
