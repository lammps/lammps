
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# This is a tutorial, not a test, so we don't ask CTest to run it.

IF(Kokkos_ENABLE_CXX11)

TRIBITS_ADD_EXECUTABLE(
  tutorial_hierarchicalparallelism_03_vectorization
  SOURCES vectorization.cpp
  COMM serial mpi
  )

ENDIF()

