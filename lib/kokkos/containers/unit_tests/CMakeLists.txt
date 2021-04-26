
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src )
KOKKOS_INCLUDE_DIRECTORIES(${KOKKOS_SOURCE_DIR}/core/unit_test/category_files)

foreach(Tag Threads;Serial;OpenMP;HPX;Cuda;HIP;SYCL)
  # Because there is always an exception to the rule
  if(Tag STREQUAL "Threads")
    set(DEVICE "PTHREAD")
  else()
    string(TOUPPER ${Tag} DEVICE)
  endif()
  string(TOLOWER ${Tag} dir)
  # Add test for that backend if it is enabled
  if(Kokkos_ENABLE_${DEVICE})
    set(UnitTestSources UnitTestMain.cpp)
    set(dir ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    file(MAKE_DIRECTORY ${dir})
    foreach(Name
        Bitset
        DualView
        DynamicView
        DynViewAPI_generic
        DynViewAPI_rank12345
        DynViewAPI_rank67
        ErrorReporter
        OffsetView
        ScatterView
        StaticCrsGraph
        UnorderedMap
        Vector
        ViewCtorPropEmbeddedDim
        )
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      file(WRITE ${dir}/dummy.cpp
          "#include <Test${Tag}_Category.hpp>\n"
          "#include <Test${Name}.hpp>\n"
      )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND UnitTestSources ${file})
    endforeach()
    KOKKOS_ADD_EXECUTABLE_AND_TEST(UnitTest_${Tag} SOURCES ${UnitTestSources})
  endif()
endforeach()
