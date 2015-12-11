
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (Kokkos_ENABLE_CXX11)
  # This is a tutorial, not a test, so we don't ask CTest to run it.
  TRIBITS_ADD_EXECUTABLE(
    tutorial_03_simple_view_lambda
    SOURCES simple_view_lambda.cpp
    COMM serial mpi
    )
ENDIF ()
