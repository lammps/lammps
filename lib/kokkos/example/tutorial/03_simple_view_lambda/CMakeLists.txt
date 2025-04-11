kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkos_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# This is a tutorial, not a test, so we don't ask CTest to run it.
kokkos_add_executable(tutorial_03_simple_view_lambda SOURCES simple_view_lambda.cpp)
