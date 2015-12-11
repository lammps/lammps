
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# This is a tutorial, not a test, so we don't ask CTest to run it.
TRIBITS_ADD_EXECUTABLE(
  tutorial_advancedviews_02_memory_traits
  SOURCES memory_traits.cpp
  COMM serial mpi
  )
