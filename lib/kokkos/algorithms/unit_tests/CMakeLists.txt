
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

    # -------------------------
    # Sort1d,3d, Random
    # -------------------------
    set(SOURCES_A)
    if(Tag STREQUAL "OpenMP")
      LIST(APPEND SOURCES_A
	TestOpenMP_Sort1D.cpp
	TestOpenMP_Sort3D.cpp
	TestOpenMP_SortDynamicView.cpp
	)
    endif()

    # Each of these inputs is an .hpp file.
    # Generate a .cpp file for each one that runs it on the current backend (Tag),
    # and add this .cpp file to the sources for UnitTest_RandomAndSort.
    foreach(SOURCES_A_Input
        TestRandomCommon
        TestSortCommon
        TestNestedSort
    )
      set(file ${dir}/${SOURCES_A_Input}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp
        "#include <Test${Tag}_Category.hpp>\n"
        "#include <${SOURCES_A_Input}.hpp>\n"
        )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND SOURCES_A ${file})
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

# FIXME_OPENMPTARGET This test causes internal compiler errors as of 09/01/22
# when compiling for Intel's Xe-HP GPUs.
if(NOT (KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL IntelLLVM))
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_RandomAndSort
    SOURCES
      UnitTestMain.cpp
      ${SOURCES_A}
  )
endif()

foreach(ID A;B;C;D;E)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_StdSet_${ID}
    SOURCES
    UnitTestMain.cpp
    ${STDALGO_SOURCES_${ID}}
    )
endforeach()

# FIXME_OPENMPTARGET This test causes internal compiler errors as of 09/01/22
# when compiling for Intel's Xe-HP GPUs.
if(NOT (KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL IntelLLVM))
  KOKKOS_ADD_EXECUTABLE(
    UnitTest_StdAlgoCompileOnly
    SOURCES TestStdAlgorithmsCompileOnly.cpp
  )
endif()
