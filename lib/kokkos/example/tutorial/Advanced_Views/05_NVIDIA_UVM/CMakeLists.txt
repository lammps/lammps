
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (Kokkos_ENABLE_Cuda_UVM)
  # This is a tutorial, not a test, so we don't ask CTest to run it.
  TRIBITS_ADD_EXECUTABLE(
    tutorial_advancedviews_05_nvidia_uvm
    SOURCES uvm_example.cpp
    COMM serial mpi
    DEPLIBS kokkoscontainers kokkoscore
    )
ENDIF ()
