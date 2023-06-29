
#Leave these here for now - I don't need transitive deps anyway
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src )
KOKKOS_INCLUDE_DIRECTORIES(${KOKKOS_SOURCE_DIR}/core/unit_test/category_files)


SET(ALGORITHM UnitTestMain.cpp)

foreach(Tag Threads;Serial;OpenMP;Cuda;HPX;HIP;SYCL;OpenMPTarget)
  string(TOUPPER ${Tag} DEVICE)
  string(TOLOWER ${Tag} dir)

  if(Kokkos_ENABLE_${DEVICE})
    set(dir ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    file(MAKE_DIRECTORY ${dir})

    # ------------------------------------------
    # Sort
    # ------------------------------------------
    # Each of these inputs is an .hpp file.
    # Generate a .cpp file for each one that runs it on the current backend (Tag),
    # and add this .cpp file to the sources for UnitTest_RandomAndSort.
    set(ALGO_SORT_SOURCES)
    foreach(SOURCE_Input
	TestSort
	TestBinSortA
	TestBinSortB
	TestNestedSort
      )
      set(file ${dir}/${SOURCE_Input}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp
        "#include <Test${Tag}_Category.hpp>\n"
        "#include <${SOURCE_Input}.hpp>\n"
        )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ALGO_SORT_SOURCES ${file})
    endforeach()

    # ------------------------------------------
    # Random
    # ------------------------------------------
    # do as above
    set(ALGO_RANDOM_SOURCES)
    foreach(SOURCE_Input
	TestRandom
      )
      set(file ${dir}/${SOURCE_Input}.cpp)
      file(WRITE ${dir}/dummy.cpp
        "#include <Test${Tag}_Category.hpp>\n"
        "#include <${SOURCE_Input}.hpp>\n"
        )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ALGO_RANDOM_SOURCES ${file})
    endforeach()

    # ------------------------------------------
    # std set A
    # ------------------------------------------
    set(STDALGO_SOURCES_A)
    foreach(Name
	StdReducers
	StdAlgorithmsConstraints
	RandomAccessIterator
	)
      list(APPEND STDALGO_SOURCES_A Test${Name}.cpp)
    endforeach()

    # ------------------------------------------
    # std set B
    # ------------------------------------------
    set(STDALGO_SOURCES_B)
    foreach(Name
	StdAlgorithmsCommon
	StdAlgorithmsMinMaxElementOps
	)
      list(APPEND STDALGO_SOURCES_B Test${Name}.cpp)
    endforeach()

    # ------------------------------------------
    # std set C
    # ------------------------------------------
    set(STDALGO_SOURCES_C)
    foreach(Name
	StdAlgorithmsCommon
	StdAlgorithmsLexicographicalCompare
	StdAlgorithmsForEach
	StdAlgorithmsFind
	StdAlgorithmsFindFirstOf
	StdAlgorithmsFindEnd
	StdAlgorithmsCount
	StdAlgorithmsEqual
	StdAlgorithmsAllAnyNoneOf
	StdAlgorithmsAdjacentFind
	StdAlgorithmsSearch
	StdAlgorithmsSearch_n
	StdAlgorithmsMismatch
	StdAlgorithmsMoveBackward
	)
      list(APPEND STDALGO_SOURCES_C Test${Name}.cpp)
    endforeach()

    # ------------------------------------------
    # std set D
    # ------------------------------------------
    set(STDALGO_SOURCES_D)
    foreach(Name
	StdAlgorithmsCommon
	StdAlgorithmsModOps
	StdAlgorithmsModSeqOps
	StdAlgorithmsReplace
	StdAlgorithmsReplaceIf
	StdAlgorithmsReplaceCopy
	StdAlgorithmsReplaceCopyIf
	StdAlgorithmsCopyIf
	StdAlgorithmsUnique
	StdAlgorithmsUniqueCopy
	StdAlgorithmsRemove
	StdAlgorithmsRemoveIf
	StdAlgorithmsRemoveCopy
	StdAlgorithmsRemoveCopyIf
	StdAlgorithmsRotate
	StdAlgorithmsRotateCopy
	StdAlgorithmsReverse
	StdAlgorithmsShiftLeft
	StdAlgorithmsShiftRight
	)
      list(APPEND STDALGO_SOURCES_D Test${Name}.cpp)
    endforeach()

    # ------------------------------------------
    # std set E
    # ------------------------------------------
    set(STDALGO_SOURCES_E)
    foreach(Name
	StdAlgorithmsCommon
	StdAlgorithmsIsSorted
	StdAlgorithmsIsSortedUntil
	StdAlgorithmsPartitioningOps
	StdAlgorithmsPartitionCopy
	StdAlgorithmsNumerics
	StdAlgorithmsAdjacentDifference
	StdAlgorithmsExclusiveScan
	StdAlgorithmsInclusiveScan
	StdAlgorithmsTransformUnaryOp
	StdAlgorithmsTransformExclusiveScan
	StdAlgorithmsTransformInclusiveScan
	)
      list(APPEND STDALGO_SOURCES_E Test${Name}.cpp)
    endforeach()

  endif()
endforeach()

# FIXME_OPENMPTARGET This test causes internal compiler errors as of 09/01/22
# when compiling for Intel's Xe-HP GPUs.
# FRIZZI: 04/26/2023: not sure if the compilation error is still applicable
# but we conservatively leave this guard on
if(NOT (KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL IntelLLVM))
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Sort
    SOURCES
    UnitTestMain.cpp
    ${ALGO_SORT_SOURCES}
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Random
    SOURCES
    UnitTestMain.cpp
    ${ALGO_RANDOM_SOURCES}
  )
endif()

# FIXME_OPENMPTARGET These tests cause internal compiler errors as of 09/01/22
# when compiling for Intel's Xe-HP GPUs.
if(KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL IntelLLVM)
  list(REMOVE_ITEM STDALGO_SOURCES_D
    TestStdAlgorithmsCopyIf.cpp
    TestStdAlgorithmsRemoveCopy.cpp
    TestStdAlgorithmsUnique.cpp
    TestStdAlgorithmsUniqueCopy.cpp
  )
  list(REMOVE_ITEM STDALGO_SOURCES_E
    TestStdAlgorithmsExclusiveScan.cpp
    TestStdAlgorithmsInclusiveScan.cpp
  )
endif()

foreach(ID A;B;C;D;E)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    AlgorithmsUnitTest_StdSet_${ID}
    SOURCES
    UnitTestMain.cpp
    ${STDALGO_SOURCES_${ID}}
    )
endforeach()

# FIXME_OPENMPTARGET This test causes internal compiler errors as of 09/01/22
# when compiling for Intel's Xe-HP GPUs.
if(NOT (KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL IntelLLVM))
  KOKKOS_ADD_EXECUTABLE(
    AlgorithmsUnitTest_StdAlgoCompileOnly
    SOURCES TestStdAlgorithmsCompileOnly.cpp
  )
endif()
