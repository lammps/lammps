
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src )

foreach(Tag Threads;Serial;OpenMP;HPX;Cuda)
  # Because there is always an exception to the rule
  if(Tag STREQUAL "Threads")
    set(DEVICE "PTHREAD")
  else()
    string(TOUPPER ${Tag} DEVICE)
  endif()
  string(TOLOWER ${Tag} dir)
  # Add test for that backend if it is enabled
  if(Kokkos_ENABLE_${DEVICE})
    KOKKOS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_${Tag}
      SOURCES
        UnitTestMain.cpp
        ${dir}/Test${Tag}_BitSet.cpp
        ${dir}/Test${Tag}_DualView.cpp
        ${dir}/Test${Tag}_DynamicView.cpp
        ${dir}/Test${Tag}_DynRankViewAPI_generic.cpp
        ${dir}/Test${Tag}_DynRankViewAPI_rank12345.cpp
        ${dir}/Test${Tag}_DynRankViewAPI_rank67.cpp
        ${dir}/Test${Tag}_ErrorReporter.cpp
        ${dir}/Test${Tag}_OffsetView.cpp
        ${dir}/Test${Tag}_ScatterView.cpp
        ${dir}/Test${Tag}_StaticCrsGraph.cpp
        ${dir}/Test${Tag}_UnorderedMap.cpp
        ${dir}/Test${Tag}_Vector.cpp
        ${dir}/Test${Tag}_ViewCtorPropEmbeddedDim.cpp
      )
  endif()
endforeach()
