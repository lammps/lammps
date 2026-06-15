include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# This is a tutorial, not a test, so we don't ask CTest to run it.
kokkos_add_executable(tutorial_advancedviews_07_overlapping_deepcopy SOURCES overlapping_deepcopy.cpp)
