#
# Add test-only library for gtest to be reused by all the subpackages
#

IF(NOT GTest_FOUND)  # fallback to internal gtest
  SET(GTEST_SOURCE_DIR ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/tpls/gtest)

  #need here for tribits
  KOKKOS_INCLUDE_DIRECTORIES(${GTEST_SOURCE_DIR})
  KOKKOS_ADD_TEST_LIBRARY(
    kokkos_gtest
    HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
    SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
  )

  TARGET_INCLUDE_DIRECTORIES(kokkos_gtest PUBLIC ${GTEST_SOURCE_DIR})
  IF((NOT (Kokkos_ENABLE_CUDA AND WIN32)) AND (NOT ("${KOKKOS_CXX_COMPILER_ID}" STREQUAL "Fujitsu")))
    TARGET_COMPILE_FEATURES(kokkos_gtest PUBLIC cxx_std_14)
  ENDIF()

  # Suppress clang-tidy diagnostics on code that we do not have control over
  IF(CMAKE_CXX_CLANG_TIDY)
    SET_TARGET_PROPERTIES(kokkos_gtest PROPERTIES CXX_CLANG_TIDY "")
  ENDIF()

  FIND_PACKAGE(Threads QUIET)
  IF(TARGET Threads::Threads)
    SET_TARGET_PROPERTIES(kokkos_gtest PROPERTIES
                          INTERFACE_LINK_LIBRARIES Threads::Threads)
  ENDIF()
ENDIF()

#
# Define Incremental Testing Feature Levels
# Define Device name mappings (i.e. what comes after Kokkos:: for the ExecSpace)
#

SET(KOKKOS_CUDA_FEATURE_LEVEL 999)
SET(KOKKOS_CUDA_NAME Cuda)
SET(KOKKOS_HIP_FEATURE_LEVEL 999)
SET(KOKKOS_HIP_NAME Experimental::HIP)
SET(KOKKOS_HPX_FEATURE_LEVEL 999)
SET(KOKKOS_HPX_NAME Experimental::HPX)
SET(KOKKOS_OPENMP_FEATURE_LEVEL 999)
SET(KOKKOS_OPENMP_NAME OpenMP)

# FIXME_OPENMPTARGET - The NVIDIA HPC compiler nvc++ only compiles the first 8 incremental tests for the OpenMPTarget backend.
IF(KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  SET(KOKKOS_OPENMPTARGET_FEATURE_LEVEL 10)
ELSE()
  SET(KOKKOS_OPENMPTARGET_FEATURE_LEVEL 14)
ENDIF()

SET(KOKKOS_OPENMPTARGET_NAME Experimental::OpenMPTarget)
SET(KOKKOS_SERIAL_FEATURE_LEVEL 999)
SET(KOKKOS_SERIAL_NAME Serial)
SET(KOKKOS_SYCL_FEATURE_LEVEL 999)
SET(KOKKOS_SYCL_NAME Experimental::SYCL)
SET(KOKKOS_THREADS_FEATURE_LEVEL 999)
SET(KOKKOS_THREADS_NAME Threads)
SET(KOKKOS_OPENACC_FEATURE_LEVEL 3)
SET(KOKKOS_OPENACC_NAME Experimental::OpenACC)


#
# Define the tests
#

#I will leave these alone for now because I don't need transitive dependencies on tests
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
KOKKOS_INCLUDE_DIRECTORIES(${KOKKOS_SOURCE_DIR}/core/unit_test/category_files)

SET(COMPILE_ONLY_SOURCES
  TestArray.cpp
  TestCreateMirror.cpp
  TestDetectionIdiom.cpp
  TestInterOp.cpp
  TestLegionInteroperability.cpp
  TestStringManipulation.cpp
  TestTypeList.cpp
)
# TestInterOp has a dependency on containers
IF(KOKKOS_HAS_TRILINOS)
  LIST(REMOVE_ITEM COMPILE_ONLY_SOURCES TestInterOp.cpp)
ENDIF()
KOKKOS_ADD_EXECUTABLE(
  TestCompileOnly
  SOURCES
  TestCompileMain.cpp
  ${COMPILE_ONLY_SOURCES}
)

foreach(Tag Threads;Serial;OpenMP;Cuda;HPX;OpenMPTarget;HIP;SYCL)
  string(TOUPPER ${Tag} DEVICE)
  string(TOLOWER ${Tag} dir)

  if(Kokkos_ENABLE_${DEVICE})
    set(dir ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    file(MAKE_DIRECTORY ${dir})
    # Needed to split this for Windows NVCC, since it ends up putting everything on the
    # command line in an intermediate compilation step even if CMake generated a response
    # file. That then exceeded the shell command line max length.
    set(${Tag}_SOURCES1A)
    foreach(Name
        AtomicOperations_int
        AtomicOperations_unsignedint
        AtomicOperations_longint
        AtomicOperations_unsignedlongint
        AtomicOperations_longlongint
        AtomicOperations_double
        AtomicOperations_float
        AtomicOperations_complexdouble
        AtomicOperations_complexfloat
        AtomicOperations_shared
        AtomicViews
        Atomics
        BlockSizeDeduction
        Concepts
        Complex
        Crs
        DeepCopyAlignment
        ExecutionSpace
        FunctorAnalysis
        Init
        JoinBackwardCompatibility
        LocalDeepCopy
        MinMaxClamp
        MathematicalConstants
        MathematicalFunctions1
        MathematicalFunctions2
        MDRange_a
        MDRange_b
        MDRange_c
        HostSharedPtr
        HostSharedPtrAccessOnDevice
        QuadPrecisionMath
        ExecSpacePartitioning
        MathematicalSpecialFunctions
        )
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp
          "#include <Test${Tag}_Category.hpp>\n"
          "#include <Test${Name}.hpp>\n"
      )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES1A ${file})
    endforeach()

    set(${Tag}_SOURCES1B)
    foreach(Name
        MDRange_d
        MDRange_e
        MDRange_f
        MDRange_g
        NumericTraits
        Other
        RangePolicy
        RangePolicyRequire
        Reductions
        Reducers_a
        Reducers_b
        Reducers_c
        Reducers_d
        Reducers_e
        Reductions_DeviceView
        Scan
        SharedAlloc
        ViewMapping_a
        )
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp
          "#include <Test${Tag}_Category.hpp>\n"
          "#include <Test${Name}.hpp>\n"
      )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES1B ${file})
    endforeach()

    SET(SOURCES2A_NAME_LIST
      Abort
      TeamBasic
      TeamReductionScan
      TeamScan
      TeamScratch
      TeamTeamSize
      TeamVectorRange
      UniqueToken
      ViewAPI_a
      ViewAPI_b
      ViewAPI_c
      ViewAPI_d
      ViewAPI_e
      ViewCopy_a
      ViewCopy_b
      ViewCtorDimMatch
      ViewHooks
      ViewLayoutStrideAssignment
      ViewMapping_b
      ViewMapping_subview
      ViewMemoryAccessViolation
      ViewOfClass
      ViewResize
      View_64bit
      WorkGraph
      WithoutInitializing
    )
    IF(KOKKOS_HAS_TRILINOS)
      LIST(REMOVE_ITEM SOURCES2A_NAME_LIST
        Abort
        ViewMemoryAccessViolation
      )
    ENDIF()

    SET(${Tag}_SOURCES2A)
    foreach(Name ${SOURCES2A_NAME_LIST})
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp
          "#include <Test${Tag}_Category.hpp>\n"
          "#include <Test${Name}.hpp>\n"
      )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES2A ${file})
    endforeach()

    set(TagHostAccessible ${Tag})
    if (Tag STREQUAL "Cuda")
      set(TagHostAccessible CudaUVM)
    elseif(Tag STREQUAL "HIP")
      set(TagHostAccessible HIPManaged)
    elseif(Tag STREQUAL "SYCL")
      set(TagHostAccessible SYCLSharedUSM)
    endif()

    set(${Tag}_SOURCES2B)
    foreach(Name
      SubView_a
      SubView_b
      SubView_c01
      SubView_c02
      SubView_c03
      SubView_c04
      SubView_c05
      )
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp
          "#include <Test${TagHostAccessible}_Category.hpp>\n"
          "#include <Test${Name}.hpp>\n"
      )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES2B ${file})
    endforeach()

    set(${Tag}_SOURCES2C)
    foreach(Name
      SubView_c06
      SubView_c07
      SubView_c08
      SubView_c09
      )
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp
          "#include <Test${TagHostAccessible}_Category.hpp>\n"
          "#include <Test${Name}.hpp>\n"
      )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES2C ${file})
    endforeach()

    set(${Tag}_SOURCES2D)
    foreach(Name
      SubView_c10
      SubView_c11
      SubView_c12
      SubView_c13
      SubView_c14
      )
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp
          "#include <Test${TagHostAccessible}_Category.hpp>\n"
          "#include <Test${Name}.hpp>\n"
      )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES2D ${file})
    endforeach()

    SET(${Tag}_SOURCES1 ${${Tag}_SOURCES1A} ${${Tag}_SOURCES1B})
    SET(${Tag}_SOURCES2 ${${Tag}_SOURCES2A} ${${Tag}_SOURCES2B} ${${Tag}_SOURCES2C} ${${Tag}_SOURCES2D})
    SET(${Tag}_SOURCES ${${Tag}_SOURCES1} ${${Tag}_SOURCES2})
  endif()
endforeach()

foreach(PairDeviceSpace HIP-HostPinned;HIP-Managed;Cuda-HostPinned;Cuda-UVM;SYCL-HostUSM;SYCL-SharedUSM)
  string(REGEX REPLACE "([^-]*)-(.*)" "\\1" DEVICE ${PairDeviceSpace})
  string(REGEX REPLACE "([^-]*)-(.*)" "\\2" SPACE ${PairDeviceSpace})

  string(TOUPPER ${DEVICE} UPPER_DEVICE)
  string(TOLOWER ${DEVICE} dir)

  if(Kokkos_ENABLE_${UPPER_DEVICE})
    set(dir ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    file(MAKE_DIRECTORY ${dir})
    foreach(Name
      SharedAlloc
      ViewAPI_a
      ViewAPI_b
      ViewAPI_c
      ViewAPI_d
      ViewAPI_e
      ViewCopy_a
      ViewCopy_b
      ViewMapping_a
      ViewMapping_b
      ViewMapping_subview
      )
      set(file ${dir}/Test${DEVICE}${SPACE}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp
          "#include <Test${DEVICE}${SPACE}_Category.hpp>\n"
          "#include <Test${Name}.hpp>\n"
      )
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${DEVICE}_SOURCES3 ${file})
    endforeach()
    list(APPEND ${DEVICE}_SOURCES ${${DEVICE}_SOURCES3})
  endif()
endforeach()

if(Kokkos_ENABLE_OPENMPTARGET)
  list(REMOVE_ITEM OpenMPTarget_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_complexfloat.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_complexdouble.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Crs.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_LocalDeepCopy.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Other.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamReductionScan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamScan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_e.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewCopy_a.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewCopy_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewMapping_subview.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewOfClass.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_WorkGraph.cpp
    )
endif()

# FIXME_OPENMPTARGET - Comment non-passing tests with amdclang++
IF(KOKKOS_ARCH_VEGA906 OR KOKKOS_ARCH_VEGA908 OR KOKKOS_ARCH_VEGA90A)
  SET(KOKKOS_AMDGPU_ARCH TRUE)
ENDIF()
IF(KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL Clang AND KOKKOS_AMDGPU_ARCH)
  LIST(REMOVE_ITEM OpenMPTarget_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_e.cpp
  )
ENDIF()

# FIXME_OPENMPTARGET - Comment non-passing tests with the NVIDIA HPC compiler nvc++
IF(KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  list(REMOVE_ITEM OpenMPTarget_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_int64_t_reduce.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_int64_t_reduce_dynamic.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_int64_t_reduce_dynamic_view.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_double_reduce.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_double_reduce_dynamic.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamTeamSize.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reductions_DeviceView.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamVectorRange.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_UniqueToken.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_HostSharedPtr.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_HostSharedPtrAccessOnDevice.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamScratch.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TestScan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TestTeamScan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TestTeamReductionScan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Atomics.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_float.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_int.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_longint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_longlongint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_double.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_unsignedint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_unsignedlongint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicViews.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_BlockSizeDeduction.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_a.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_c.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_d.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_e.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewMapping_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamBasic.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Scan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_NumericTraits.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_DeepCopyAlignment.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MathematicalFunctions.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c01.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c02.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c03.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c04.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c05.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c06.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c07.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c08.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c09.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c10.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c11.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c12.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c13.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MDRange_a.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MDRange_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MDRange_c.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MDRange_d.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_a.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_c.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_d.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_f.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewResize.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_RangePolicyRequire.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_RangePolicy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/default/TestDefaultDeviceType_a1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/default/TestDefaultDeviceType_b1.cpp
    )
endif()

if(Kokkos_ENABLE_SERIAL)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Serial1
    SOURCES
    UnitTestMainInit.cpp
    ${Serial_SOURCES1}
    serial/TestSerial_Task.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Serial2
    SOURCES
    UnitTestMainInit.cpp
    ${Serial_SOURCES2}
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SerialGraph
    SOURCES
    UnitTestMainInit.cpp
    serial/TestSerial_Graph.cpp
  )
endif()

if(Kokkos_ENABLE_THREADS)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Threads
    SOURCES ${Threads_SOURCES}
    UnitTestMainInit.cpp
  )
endif()

if (Kokkos_ENABLE_OPENMP)
  set(OpenMP_EXTRA_SOURCES
    openmp/TestOpenMP_Task.cpp
    openmp/TestOpenMP_PartitionMaster.cpp
  )
  if (Kokkos_ENABLE_DEPRECATED_CODE_3)
    list(APPEND OpenMP_EXTRA_SOURCES openmp/TestOpenMP_Task.cpp)
  endif ()
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_OpenMP
    SOURCES
    UnitTestMainInit.cpp
    ${OpenMP_SOURCES}
    ${OpenMP_EXTRA_SOURCES}
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_OpenMPInterOp
    SOURCES
      UnitTestMain.cpp
      openmp/TestOpenMP_InterOp.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_OpenMPGraph
    SOURCES
      UnitTestMainInit.cpp
      openmp/TestOpenMP_Graph.cpp
  )
endif()

if(Kokkos_ENABLE_HPX)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HPX
    SOURCES
      UnitTestMainInit.cpp
      ${HPX_SOURCES}
      hpx/TestHPX_Task.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HPXInterOp
    SOURCES
      UnitTestMain.cpp
      hpx/TestHPX_InterOp.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HPX_IndependentInstances
    SOURCES
      UnitTestMainInit.cpp
      hpx/TestHPX_IndependentInstances.cpp
      hpx/TestHPX_IndependentInstancesDelayedExecution.cpp
      hpx/TestHPX_IndependentInstancesInstanceIds.cpp
      hpx/TestHPX_IndependentInstancesRefCounting.cpp
  )
endif()

if(Kokkos_ENABLE_OPENMPTARGET)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_OpenMPTarget
    SOURCES
    UnitTestMainInit.cpp
    ${OpenMPTarget_SOURCES}
  )
endif()

if(Kokkos_ENABLE_CUDA)
    KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Cuda1
    SOURCES
      UnitTestMainInit.cpp
      ${Cuda_SOURCES1}
      cuda/TestCuda_ReducerViewSizeLimit.cpp
    )

    KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Cuda2
    SOURCES
      UnitTestMainInit.cpp
      ${Cuda_SOURCES2}
    )

    KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Cuda3
    SOURCES
      UnitTestMainInit.cpp
      cuda/TestCuda_Task.cpp
      cuda/TestCuda_TeamScratchStreams.cpp
      ${Cuda_SOURCES3}
      cuda/TestCuda_Spaces.cpp
  )

    KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_CudaTimingBased
    SOURCES
      UnitTestMainInit.cpp
      cuda/TestCuda_DebugSerialExecution.cpp
      cuda/TestCuda_DebugPinUVMSpace.cpp
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_CudaInterOpInit
    SOURCES
      UnitTestMain.cpp
      cuda/TestCuda_InterOp_Init.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_CudaInterOpStreams
    SOURCES
      UnitTestMain.cpp
      cuda/TestCuda_InterOp_Streams.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_CudaGraph
    SOURCES
      UnitTestMainInit.cpp
      cuda/TestCuda_Graph.cpp
  )
endif()

if(Kokkos_ENABLE_HIP)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HIP
    SOURCES
      UnitTestMainInit.cpp
      ${HIP_SOURCES}
      hip/TestHIP_ScanUnit.cpp
      hip/TestHIP_Spaces.cpp
      hip/TestHIP_Memory_Requirements.cpp
      hip/TestHIP_TeamScratchStreams.cpp
      hip/TestHIP_AsyncLauncher.cpp
      hip/TestHIP_BlocksizeDeduction.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HIPInterOpInit
    SOURCES
      UnitTestMain.cpp
      hip/TestHIP_InterOp_Init.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HIPInterOpStreams
    SOURCES
      UnitTestMain.cpp
      hip/TestHIP_InterOp_Streams.cpp
  )
endif()

if(Kokkos_ENABLE_SYCL)
  list(REMOVE_ITEM SYCL_SOURCES1A
       # FIXME_SYCL atomic_fetch_oper for large types to be implemented
       ${CMAKE_CURRENT_BINARY_DIR}/sycl/TestSYCL_AtomicOperations_complexdouble.cpp
  )

  list(REMOVE_ITEM SYCL_SOURCES2A
       ${CMAKE_CURRENT_BINARY_DIR}/sycl/TestSYCL_WorkGraph.cpp
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCL1A
    SOURCES
      UnitTestMainInit.cpp
      ${SYCL_SOURCES1A}
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCL1B
    SOURCES
      UnitTestMainInit.cpp
      ${SYCL_SOURCES1B}
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCL2A
    SOURCES
      UnitTestMainInit.cpp
      ${SYCL_SOURCES2A}
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCL2B
    SOURCES
      UnitTestMainInit.cpp
      ${SYCL_SOURCES2B}
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCL2C
    SOURCES
      UnitTestMainInit.cpp
      ${SYCL_SOURCES2C}
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCL2D
    SOURCES
      UnitTestMainInit.cpp
      ${SYCL_SOURCES2D}
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCL3
    SOURCES
      UnitTestMainInit.cpp
      # FIXME_SYCL
      sycl/TestSYCL_Task.cpp
      sycl/TestSYCL_TeamScratchStreams.cpp
      ${SYCL_SOURCES3}
      sycl/TestSYCL_Spaces.cpp
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCLInterOpInit
    SOURCES
      UnitTestMain.cpp
      sycl/TestSYCL_InterOp_Init.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCLInterOpInit_Context
    SOURCES
    UnitTestMainInit.cpp
      sycl/TestSYCL_InterOp_Init_Context.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_SYCLInterOpStreams
    SOURCES
      UnitTestMain.cpp
     sycl/TestSYCL_InterOp_Streams.cpp
  )
endif()

# FIXME_OPENMPTARGET - Comment non-passing tests with the NVIDIA HPC compiler nvc++
if (KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  SET(DEFAULT_DEVICE_SOURCES
    UnitTestMainInit.cpp
    TestInitializationSettings.cpp
    TestParseCmdLineArgsAndEnvVars.cpp
    default/TestDefaultDeviceType.cpp
  )
else()
  SET(DEFAULT_DEVICE_SOURCES
    UnitTestMainInit.cpp
    TestInitializationSettings.cpp
    TestParseCmdLineArgsAndEnvVars.cpp
    default/TestDefaultDeviceType.cpp
    default/TestDefaultDeviceType_a1.cpp
    default/TestDefaultDeviceType_b1.cpp
    default/TestDefaultDeviceType_c1.cpp
    default/TestDefaultDeviceType_a2.cpp
    default/TestDefaultDeviceType_b2.cpp
    default/TestDefaultDeviceType_c2.cpp
    default/TestDefaultDeviceType_a3.cpp
    default/TestDefaultDeviceType_b3.cpp
    default/TestDefaultDeviceType_c3.cpp
    default/TestDefaultDeviceType_d.cpp
    default/TestDefaultDeviceTypeResize.cpp
    default/TestDefaultDeviceTypeViewAPI.cpp
  )
endif()

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Default
  SOURCES ${DEFAULT_DEVICE_SOURCES}
)

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_PushFinalizeHook
  SOURCES
    UnitTest_PushFinalizeHook.cpp
)

# This test is intended for development and debugging by putting code
# into TestDefaultDeviceDevelop.cpp. By default its empty.
KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Develop
  SOURCES
    UnitTestMainInit.cpp
    default/TestDefaultDeviceDevelop.cpp
)

# This test is special, because it passes exactly when it prints the
# message "PASSED: I am the custom std::terminate handler.", AND calls
# std::terminate.  This means that we can't use
# KOKKOS_ADD_EXECUTABLE_AND_TEST.  See GitHub issue #2147.

KOKKOS_ADD_TEST_EXECUTABLE( push_finalize_hook_terminate
  SOURCES UnitTest_PushFinalizeHook_terminate.cpp
)

KOKKOS_ADD_ADVANCED_TEST( UnitTest_PushFinalizeHook_terminate
  TEST_0
    EXEC push_finalize_hook_terminate
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION
      "PASSED: I am the custom std::terminate handler."
    ALWAYS_FAIL_ON_ZERO_RETURN
)
  if(KOKKOS_ENABLE_TUNING)
    KOKKOS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_TuningBuiltins
      SOURCES
      tools/TestBuiltinTuners.cpp
    )
    KOKKOS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_TuningBasics
      SOURCES
        tools/TestTuning.cpp
    )
    KOKKOS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_CategoricalTuner
      SOURCES
      tools/TestCategoricalTuner.cpp
    )
  endif()
  if(NOT Kokkos_ENABLE_OPENMPTARGET)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_LogicalSpaces
    SOURCES
      tools/TestLogicalSpaces.cpp
  )
  endif()
  if(NOT (KOKKOS_CXX_COMPILER_ID STREQUAL Intel AND KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 18.0.0))
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_KokkosP
    SOURCES
    UnitTestMainInit.cpp
    tools/TestEventCorrectness.cpp
    tools/TestWithoutInitializing.cpp
    tools/TestProfilingSection.cpp
  )
  endif()
  if(KOKKOS_ENABLE_LIBDL)
    KOKKOS_ADD_EXECUTABLE_AND_TEST(
      UnitTest_ToolIndependence
      SOURCES
      tools/TestIndependence.cpp
    )
    TARGET_COMPILE_DEFINITIONS(
      KokkosCore_UnitTest_ToolIndependence PUBLIC
      KOKKOS_TOOLS_INDEPENDENT_BUILD
    )
    KOKKOS_ADD_TEST_LIBRARY(
      kokkosprinter-tool SHARED
      SOURCES tools/printing-tool.cpp
    )

    if((NOT (Kokkos_ENABLE_CUDA AND WIN32)) AND (NOT ("${KOKKOS_CXX_COMPILER_ID}" STREQUAL "Fujitsu")))
      TARGET_COMPILE_FEATURES(kokkosprinter-tool PUBLIC cxx_std_14)
    endif()

    KOKKOS_ADD_TEST_EXECUTABLE(
      ProfilingAllCalls
      tools/TestAllCalls.cpp
    )

    set(ADDRESS_REGEX "0x[0-9a-f]*")
    set(MEMSPACE_REGEX "[HC][ou][sd][ta][a-zA-Z]*")
    set(SIZE_REGEX "[0-9]*")
    set(SKIP_SCRATCH_INITIALIZATION_REGEX ".*")

    # check help works via environment variable
    KOKKOS_ADD_TEST(
      SKIP_TRIBITS
      NAME ProfilingTestLibraryLoadHelp
      EXE  ProfilingAllCalls
      TOOL kokkosprinter-tool
      ARGS --kokkos-tools-help
      PASS_REGULAR_EXPRESSION
        "kokkosp_init_library::kokkosp_print_help:KokkosCore_ProfilingAllCalls::kokkosp_finalize_library::")

    # check help works via direct library specification
    KOKKOS_ADD_TEST(
      SKIP_TRIBITS
      NAME ProfilingTestLibraryCmdLineHelp
      EXE  ProfilingAllCalls
      ARGS --kokkos-tools-help
           --kokkos-tools-libs=$<TARGET_FILE:kokkosprinter-tool>
      PASS_REGULAR_EXPRESSION
        "kokkosp_init_library::kokkosp_print_help:KokkosCore_ProfilingAllCalls::kokkosp_finalize_library::")

    KOKKOS_ADD_TEST(
      SKIP_TRIBITS
      NAME ProfilingTestLibraryLoad
      EXE  ProfilingAllCalls
      TOOL kokkosprinter-tool
      ARGS --kokkos-tools-args="-c test delimit"
      PASS_REGULAR_EXPRESSION "kokkosp_init_library::kokkosp_parse_args:4:KokkosCore_ProfilingAllCalls:-c:test:delimit::.*::kokkosp_allocate_data:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::kokkosp_begin_parallel_for:Kokkos::View::initialization [[]source] via memset:[0-9]+:0::kokkosp_end_parallel_for:0::kokkosp_allocate_data:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:40::kokkosp_begin_parallel_for:Kokkos::View::initialization [[]destination] via memset:[0-9]+:0::kokkosp_end_parallel_for:0::kokkosp_begin_deep_copy:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::.*kokkosp_end_deep_copy::kokkosp_begin_parallel_for:parallel_for:${SIZE_REGEX}:0::kokkosp_end_parallel_for:0::kokkosp_begin_parallel_reduce:parallel_reduce:${SIZE_REGEX}:1${SKIP_SCRATCH_INITIALIZATION_REGEX}::kokkosp_end_parallel_reduce:1::kokkosp_begin_parallel_scan:parallel_scan:${SIZE_REGEX}:2::kokkosp_end_parallel_scan:2::kokkosp_push_profile_region:push_region::kokkosp_pop_profile_region::kokkosp_create_profile_section:created_section:3::kokkosp_start_profile_section:3::kokkosp_stop_profile_section:3::kokkosp_destroy_profile_section:3::kokkosp_profile_event:profiling_event::kokkosp_declare_metadata:dogs:good::kokkosp_deallocate_data:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:40::kokkosp_deallocate_data:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::kokkosp_finalize_library::"
    )

    # Above will test that leading/trailing quotes are stripped bc ctest cmd args is:
    #       "--kokkos-tools-args="-c test delimit""
    # The bracket argument syntax: [=[ and ]=] used below ensures it is treated as
    # a single argument:
    #       "--kokkos-tools-args=-c test delimit"
    #
    # https://cmake.org/cmake/help/latest/manual/cmake-language.7.html#bracket-argument
    #
    KOKKOS_ADD_TEST(
      SKIP_TRIBITS
      NAME ProfilingTestLibraryCmdLine
      EXE  ProfilingAllCalls
      ARGS [=[--kokkos-tools-args=-c test delimit]=]
            --kokkos-tools-libs=$<TARGET_FILE:kokkosprinter-tool>
      PASS_REGULAR_EXPRESSION "kokkosp_init_library::kokkosp_parse_args:4:KokkosCore_ProfilingAllCalls:-c:test:delimit::.*::kokkosp_allocate_data:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::kokkosp_begin_parallel_for:Kokkos::View::initialization [[]source] via memset:[0-9]+:0::kokkosp_end_parallel_for:0::kokkosp_allocate_data:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:40::kokkosp_begin_parallel_for:Kokkos::View::initialization [[]destination] via memset:[0-9]+:0::kokkosp_end_parallel_for:0::kokkosp_begin_deep_copy:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::.*kokkosp_end_deep_copy::kokkosp_begin_parallel_for:parallel_for:${SIZE_REGEX}:0::kokkosp_end_parallel_for:0::kokkosp_begin_parallel_reduce:parallel_reduce:${SIZE_REGEX}:1${SKIP_SCRATCH_INITIALIZATION_REGEX}::kokkosp_end_parallel_reduce:1::kokkosp_begin_parallel_scan:parallel_scan:${SIZE_REGEX}:2::kokkosp_end_parallel_scan:2::kokkosp_push_profile_region:push_region::kokkosp_pop_profile_region::kokkosp_create_profile_section:created_section:3::kokkosp_start_profile_section:3::kokkosp_stop_profile_section:3::kokkosp_destroy_profile_section:3::kokkosp_profile_event:profiling_event::kokkosp_declare_metadata:dogs:good::kokkosp_deallocate_data:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:40::kokkosp_deallocate_data:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::kokkosp_finalize_library::"
    )
  endif() #KOKKOS_ENABLE_LIBDL
if(NOT KOKKOS_HAS_TRILINOS)
KOKKOS_ADD_TEST_EXECUTABLE(
  StackTraceTestExec
  SOURCES
    TestStackTrace.cpp
    TestStackTrace_f0.cpp
    TestStackTrace_f1.cpp
    TestStackTrace_f2.cpp
    TestStackTrace_f3.cpp
    TestStackTrace_f4.cpp
)
# We need -rdynamic on GNU platforms for the stacktrace functionality
# to work correctly with shared libraries
KOKKOS_SET_EXE_PROPERTY(StackTraceTestExec ENABLE_EXPORTS ON)

KOKKOS_ADD_TEST( NAME UnitTest_StackTraceTest
                 EXE  StackTraceTestExec
                 FAIL_REGULAR_EXPRESSION "FAILED"
               )
endif()

if(Kokkos_ENABLE_DEPRECATED_CODE_3)
  foreach(INITTESTS_NUM RANGE 1 18)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_DefaultInit_${INITTESTS_NUM}
    SOURCES UnitTestMain.cpp default/TestDefaultDeviceTypeInit_${INITTESTS_NUM}.cpp
  )
  endforeach(INITTESTS_NUM)
endif()

if (KOKKOS_ENABLE_HWLOC)
KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_HWLOC
  SOURCES UnitTestMain.cpp  TestHWLOC.cpp
)
endif()

FUNCTION (KOKKOS_ADD_INCREMENTAL_TEST DEVICE)
  KOKKOS_OPTION( ${DEVICE}_EXCLUDE_TESTS "" STRING "Incremental test exclude list" )
  # Add unit test main
  SET(${DEVICE}_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMainInit.cpp)

  # Iterate over incremental tests in directory

  APPEND_GLOB(INCREMENTAL_FILE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/incremental/*.hpp)

  SET(DEVICE_NAME ${KOKKOS_${DEVICE}_NAME})
  FOREACH (CURRENT_FILE_PATH ${INCREMENTAL_FILE_LIST})
    GET_FILENAME_COMPONENT( CURRENT_FILE_NAME ${CURRENT_FILE_PATH} NAME )
    STRING (REPLACE ".hpp" "" CURRENT_TEST_NAME ${CURRENT_FILE_NAME})
    IF (NOT CURRENT_TEST_NAME IN_LIST Kokkos_${DEVICE}_EXCLUDE_TESTS)
       SET (CURRENT_TEST_OUTPUT_FILENAME ${CURRENT_TEST_NAME}_${DEVICE})
       FILE( STRINGS ${CURRENT_FILE_PATH} CURRENT_REQUIRED_FEATURE_LINE REGEX "Kokkos_Feature_Level_Required" )
       # From each test get level implementation required
       STRING( REGEX REPLACE ".*Kokkos_Feature_Level_Required:" "" CURRENT_REQUIRED_FEATURE_LEVEL ${CURRENT_REQUIRED_FEATURE_LINE} )
       # Cross-reference list of dependencies with selected feature list > matching feature test files are added to test applications
       IF (KOKKOS_${DEVICE}_FEATURE_LEVEL GREATER_EQUAL CURRENT_REQUIRED_FEATURE_LEVEL)
          CONFIGURE_FILE (IncrementalTest.cpp.in ${CMAKE_BINARY_DIR}/core/unit_test/generated/${CURRENT_TEST_OUTPUT_FILENAME}.cpp )
          SET(${DEVICE}_SOURCES ${${DEVICE}_SOURCES}; ${CMAKE_BINARY_DIR}/core/unit_test/generated/${CURRENT_TEST_OUTPUT_FILENAME}.cpp)
       ENDIF()
     ENDIF()
  ENDFOREACH()

  STRING(TOUPPER ${DEVICE} UC_DEVICE)

  KOKKOS_OPTION (
    ENABLE_${UC_DEVICE} ON BOOL "ENABLE ${UC_DEVICE}"
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    IncrementalTest_${DEVICE}
    SOURCES ${${DEVICE}_SOURCES}
  )

  SET(EXE_NAME ${PACKAGE_NAME}_IncrementalTest_${DEVICE})
  # Check that the target was actually created because in a TribITS build
  # where only tests marked as PERFORMANCE enabled it would not be.
  IF(TARGET ${EXE_NAME})
    TARGET_INCLUDE_DIRECTORIES(${EXE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/incremental )
  ENDIF()

ENDFUNCTION()

FOREACH (DEVICE ${KOKKOS_ENABLED_DEVICES})
  KOKKOS_ADD_INCREMENTAL_TEST(${DEVICE})
ENDFOREACH()

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_CTestDevice
  SOURCES UnitTestMain.cpp  TestCTestDevice.cpp
)

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_CMakePassCmdLineArgs
  SOURCES UnitTest_CMakePassCmdLineArgs.cpp
  ARGS "one 2 THREE"
)

if (KOKKOS_ENABLE_HEADER_SELF_CONTAINMENT_TESTS AND NOT KOKKOS_HAS_TRILINOS AND NOT WIN32)
  add_subdirectory(headers_self_contained)
endif()
