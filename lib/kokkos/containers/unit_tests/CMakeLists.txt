kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkos_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
kokkos_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)
kokkos_include_directories(${KOKKOS_SOURCE_DIR}/core/unit_test/category_files)

foreach(Tag Threads;Serial;OpenMP;HPX;Cuda;HIP;SYCL)
  string(TOUPPER ${Tag} DEVICE)
  string(TOLOWER ${Tag} dir)
  # Add test for that backend if it is enabled
  if(KOKKOS_ENABLE_${DEVICE})
    set(UnitTestSources UnitTestMain.cpp)
    set(dir ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    file(MAKE_DIRECTORY ${dir})
    set(DeprecatedTests Vector StaticCrsGraph)
    foreach(
      Name
      Bitset
      DualView
      DynamicView
      DynViewAPI_generic
      DynViewAPI_rank12345
      DynViewAPI_rank67
      DynRankView_Ctors
      DynRankView_LayoutMember
      DynRankView_TeamScratch
      DynRankView_ViewCustomization
      ErrorReporter
      OffsetView
      ScatterView
      StaticCrsGraph
      WithoutInitializing
      UnorderedMap
      Vector
      ViewCtorPropEmbeddedDim
    )
      if(NOT Kokkos_ENABLE_DEPRECATED_CODE_4 AND Name IN_LIST DeprecatedTests)
        continue() # skip tests for deprecated features if deprecated code 4 is not enabled
      endif()
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      file(WRITE ${dir}/dummy.cpp "#include <Test${Tag}_Category.hpp>\n" "#include <Test${Name}.hpp>\n")
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND UnitTestSources ${file})
    endforeach()
    #fatal error C1128: number of sections exceeded object file format limit: compile with /bigobj
    if(KOKKOS_ENABLE_CUDA AND WIN32)
      list(REMOVE_ITEM UnitTestSources ${dir}/TestCuda_DynViewAPI_generic.cpp)
    endif()

    # FIXME_NVHPC: NVC++-S-0000-Internal compiler error. extractor: bad opc       0
    if(KOKKOS_ENABLE_CUDA AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
      list(REMOVE_ITEM UnitTestSources ${dir}/TestCuda_WithoutInitializing.cpp)
    endif()

    if(NOT Kokkos_ENABLE_IMPL_VIEW_LEGACY)
      list(REMOVE_ITEM UnitTestSources ${dir}/Test${Tag}_ViewCtorPropEmbeddedDim.cpp)
    else()
      list(REMOVE_ITEM UnitTestSources ${dir}/Test${Tag}_DynRankView_ViewCustomization.cpp)
    endif()
    kokkos_add_executable_and_test(ContainersUnitTest_${Tag} SOURCES ${UnitTestSources})
  endif()
endforeach()

set(COMPILE_ONLY_SOURCES TestCreateMirror.cpp TestDualViewParameterPack.cpp TestIsViewTrait.cpp
                         TestDynRankViewTypedefs.cpp
)
kokkos_add_executable(ContainersTestCompileOnly SOURCES TestCompileMain.cpp ${COMPILE_ONLY_SOURCES})
