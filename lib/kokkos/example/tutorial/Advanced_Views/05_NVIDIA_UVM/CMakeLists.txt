include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(Kokkos_ENABLE_CUDA_UVM)
  # This is a tutorial, not a test, so we don't ask CTest to run it.
  kokkos_add_executable(tutorial_advancedviews_05_nvidia_uvm SOURCES uvm_example.cpp)
endif()
