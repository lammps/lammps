
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# This is a tutorial, not a test, so we don't ask CTest to run it.
KOKKOS_ADD_EXECUTABLE(
  tutorial_06_simple_mdrangepolicy
  SOURCES simple_mdrangepolicy.cpp
)

