#
# Add test-only library for gtest to be reused by all the subpackages
#


SET(GTEST_SOURCE_DIR ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/tpls/gtest)

#need here for tribits
KOKKOS_INCLUDE_DIRECTORIES(${GTEST_SOURCE_DIR})
KOKKOS_ADD_TEST_LIBRARY(
  kokkos_gtest
  HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
  SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
)
#These can be direct, no need for Tribits or Kokkos wrappers

# WORKAROUND FOR HIPCC
IF(Kokkos_ENABLE_HIP)
  TARGET_COMPILE_DEFINITIONS(kokkos_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0 --amdgpu-target=gfx906")
ELSE()
  TARGET_COMPILE_DEFINITIONS(kokkos_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0")
ENDIF()

TARGET_INCLUDE_DIRECTORIES(kokkos_gtest PUBLIC ${GTEST_SOURCE_DIR})
#Gtest minimally requires C++11
TARGET_COMPILE_FEATURES(kokkos_gtest PUBLIC cxx_std_11)

#
# Define Incremental Testing Feature Levels
# Define Device name mappings (i.e. what comes after Kokkos:: for the ExecSpace)
#

SET(KOKKOS_CUDA_FEATURE_LEVEL 999)
SET(KOKKOS_CUDA_NAME Cuda)
SET(KOKKOS_HIP_FEATURE_LEVEL 12)
SET(KOKKOS_HIP_NAME Experimental::HIP)
SET(KOKKOS_HPX_FEATURE_LEVEL 999)
SET(KOKKOS_HPX_NAME Experimental::HPX)
SET(KOKKOS_OPENMP_FEATURE_LEVEL 999)
SET(KOKKOS_OPENMP_NAME OpenMP)
SET(KOKKOS_OPENMPTARGET_FEATURE_LEVEL 10)
SET(KOKKOS_OPENMPTARGET_NAME Experimental::OpenMPTarget)
SET(KOKKOS_SERIAL_FEATURE_LEVEL 999)
SET(KOKKOS_SERIAL_NAME Serial)
SET(KOKKOS_THREADS_FEATURE_LEVEL 999)
SET(KOKKOS_THREADS_NAME Threads)


#
# Define the tests
#

#I will leave these alone for now because I don't need transitive dependencies on tests
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

foreach(Tag Threads;Serial;OpenMP;Cuda;HPX;OpenMPTarget;HIP)
  # Because there is always an exception to the rule
  if(Tag STREQUAL "Threads")
    set(DEVICE "PTHREAD")
  else()
    string(TOUPPER ${Tag} DEVICE)
  endif()
  string(TOLOWER ${Tag} dir)

  SET(${Tag}_SOURCES
      UnitTestMainInit.cpp
      ${dir}/Test${Tag}_AtomicOperations_int.cpp
      ${dir}/Test${Tag}_AtomicOperations_unsignedint.cpp
      ${dir}/Test${Tag}_AtomicOperations_longint.cpp
      ${dir}/Test${Tag}_AtomicOperations_unsignedlongint.cpp
      ${dir}/Test${Tag}_AtomicOperations_longlongint.cpp
      ${dir}/Test${Tag}_AtomicOperations_double.cpp
      ${dir}/Test${Tag}_AtomicOperations_float.cpp
      ${dir}/Test${Tag}_AtomicOperations_complexdouble.cpp
      ${dir}/Test${Tag}_AtomicOperations_complexfloat.cpp
      ${dir}/Test${Tag}_AtomicViews.cpp
      ${dir}/Test${Tag}_Atomics.cpp
      ${dir}/Test${Tag}_Concepts.cpp
      ${dir}/Test${Tag}_Complex.cpp
      ${dir}/Test${Tag}_Crs.cpp
      ${dir}/Test${Tag}_DeepCopyAlignment.cpp
      ${dir}/Test${Tag}_FunctorAnalysis.cpp
      ${dir}/Test${Tag}_Init.cpp
      ${dir}/Test${Tag}_LocalDeepCopy.cpp
      ${dir}/Test${Tag}_MDRange_a.cpp
      ${dir}/Test${Tag}_MDRange_b.cpp
      ${dir}/Test${Tag}_MDRange_c.cpp
      ${dir}/Test${Tag}_MDRange_d.cpp
      ${dir}/Test${Tag}_MDRange_e.cpp
      ${dir}/Test${Tag}_Other.cpp
      ${dir}/Test${Tag}_RangePolicy.cpp
      ${dir}/Test${Tag}_RangePolicyRequire.cpp
      ${dir}/Test${Tag}_Reductions.cpp
      ${dir}/Test${Tag}_Reducers_a.cpp
      ${dir}/Test${Tag}_Reducers_b.cpp
      ${dir}/Test${Tag}_Reducers_c.cpp
      ${dir}/Test${Tag}_Reducers_d.cpp
      ${dir}/Test${Tag}_Reductions_DeviceView.cpp
      ${dir}/Test${Tag}_Scan.cpp
      ${dir}/Test${Tag}_SharedAlloc.cpp
      ${dir}/Test${Tag}_SubView_a.cpp
      ${dir}/Test${Tag}_SubView_b.cpp
      ${dir}/Test${Tag}_SubView_c01.cpp
      ${dir}/Test${Tag}_SubView_c02.cpp
      ${dir}/Test${Tag}_SubView_c03.cpp
      ${dir}/Test${Tag}_SubView_c04.cpp
      ${dir}/Test${Tag}_SubView_c05.cpp
      ${dir}/Test${Tag}_SubView_c06.cpp
      ${dir}/Test${Tag}_SubView_c07.cpp
      ${dir}/Test${Tag}_SubView_c08.cpp
      ${dir}/Test${Tag}_SubView_c09.cpp
      ${dir}/Test${Tag}_SubView_c10.cpp
      ${dir}/Test${Tag}_SubView_c11.cpp
      ${dir}/Test${Tag}_SubView_c12.cpp
      ${dir}/Test${Tag}_SubView_c13.cpp
      ${dir}/Test${Tag}_Team.cpp
      ${dir}/Test${Tag}_TeamReductionScan.cpp
      ${dir}/Test${Tag}_TeamScratch.cpp
      ${dir}/Test${Tag}_TeamTeamSize.cpp
      ${dir}/Test${Tag}_TeamVectorRange.cpp
      ${dir}/Test${Tag}_UniqueToken.cpp
      ${dir}/Test${Tag}_ViewAPI_a.cpp
      ${dir}/Test${Tag}_ViewAPI_b.cpp
      ${dir}/Test${Tag}_ViewAPI_c.cpp
      ${dir}/Test${Tag}_ViewAPI_d.cpp
      ${dir}/Test${Tag}_ViewAPI_e.cpp
      ${dir}/Test${Tag}_ViewLayoutStrideAssignment.cpp
      ${dir}/Test${Tag}_ViewMapping_a.cpp
      ${dir}/Test${Tag}_ViewMapping_b.cpp
      ${dir}/Test${Tag}_ViewMapping_subview.cpp
      ${dir}/Test${Tag}_ViewOfClass.cpp
      ${dir}/Test${Tag}_WorkGraph.cpp
      ${dir}/Test${Tag}_View_64bit.cpp
      ${dir}/Test${Tag}_ViewResize.cpp
  )
endforeach()

if(Kokkos_ENABLE_OPENMPTARGET)
  list(REMOVE_ITEM OpenMPTarget_SOURCES
    openmptarget/TestOpenMPTarget_AtomicOperations_complexdouble.cpp
    openmptarget/TestOpenMPTarget_MDRange_a.cpp
    openmptarget/TestOpenMPTarget_MDRange_b.cpp
    openmptarget/TestOpenMPTarget_MDRange_c.cpp
    openmptarget/TestOpenMPTarget_MDRange_d.cpp
    openmptarget/TestOpenMPTarget_MDRange_e.cpp
    openmptarget/TestOpenMPTarget_Other.cpp
    openmptarget/TestOpenMPTarget_Scan.cpp
    openmptarget/TestOpenMPTarget_Team.cpp
    openmptarget/TestOpenMPTarget_TeamScratch.cpp
    openmptarget/TestOpenMPTarget_ViewAPI_e.cpp
    openmptarget/TestOpenMPTarget_ViewMapping_subview.cpp
    openmptarget/TestOpenMPTarget_ViewOfClass.cpp
  )
endif()

if(Kokkos_ENABLE_HIP)
  # FIXME Linktime error: undefined reference to
  # Kokkos::Impl::ViewDimensin<0ul, ...>(unsigned int, ...)
  list(REMOVE_ITEM Serial_SOURCES serial/TestSerial_ViewLayoutStrideAssignment.cpp)
endif()

if(Kokkos_ENABLE_SERIAL)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Serial
    SOURCES
    ${Serial_SOURCES}
    serial/TestSerial_Task.cpp
  )
endif()

if(Kokkos_ENABLE_PTHREAD)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Threads
    SOURCES ${Threads_SOURCES}
  )
endif()

if(Kokkos_ENABLE_OPENMP)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_OpenMP
    SOURCES
    ${OpenMP_SOURCES}
    openmp/TestOpenMP_Task.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_OpenMPInterOp
    SOURCES
      UnitTestMain.cpp
      openmp/TestOpenMP_InterOp.cpp
  )
endif()

if(Kokkos_ENABLE_HPX)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HPX
    SOURCES
      ${HPX_SOURCES}
      hpx/TestHPX_Task.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HPXInterOp
    SOURCES
      UnitTestMain.cpp
      hpx/TestHPX_InterOp.cpp
  )
endif()

if(Kokkos_ENABLE_OPENMPTARGET)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_OpenMPTarget
    SOURCES
    ${OpenMPTarget_SOURCES}
  )
endif()

if(Kokkos_ENABLE_CUDA)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Cuda
    SOURCES
      ${Cuda_SOURCES}
      cuda/TestCuda_Task.cpp
      cuda/TestCudaHostPinned_SharedAlloc.cpp
      cuda/TestCudaHostPinned_ViewCopy.cpp
      cuda/TestCudaHostPinned_ViewAPI_a.cpp
      cuda/TestCudaHostPinned_ViewAPI_b.cpp
      cuda/TestCudaHostPinned_ViewAPI_c.cpp
      cuda/TestCudaHostPinned_ViewAPI_d.cpp
      cuda/TestCudaHostPinned_ViewAPI_e.cpp
      cuda/TestCudaHostPinned_ViewMapping_a.cpp
      cuda/TestCudaHostPinned_ViewMapping_b.cpp
      cuda/TestCudaHostPinned_ViewMapping_subview.cpp
      cuda/TestCudaUVM_SharedAlloc.cpp
      cuda/TestCudaUVM_ViewCopy.cpp
      cuda/TestCudaUVM_ViewAPI_a.cpp
      cuda/TestCudaUVM_ViewAPI_b.cpp
      cuda/TestCudaUVM_ViewAPI_c.cpp
      cuda/TestCudaUVM_ViewAPI_d.cpp
      cuda/TestCudaUVM_ViewAPI_e.cpp
      cuda/TestCudaUVM_ViewMapping_a.cpp
      cuda/TestCudaUVM_ViewMapping_b.cpp
      cuda/TestCudaUVM_ViewMapping_subview.cpp
      cuda/TestCuda_Spaces.cpp
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
endif()

if(Kokkos_ENABLE_HIP)
  # FIXME_HIP
  LIST(REMOVE_ITEM HIP_SOURCES
    hip/TestHIP_AtomicOperations_complexdouble.cpp
    hip/TestHIP_Other.cpp
    hip/TestHIP_Reductions_DeviceView.cpp
    hip/TestHIP_Team.cpp
    hip/TestHIP_TeamReductionScan.cpp
    hip/TestHIP_TeamScratch.cpp
    hip/TestHIP_TeamTeamSize.cpp
    hip/TestHIP_TeamVectorRange.cpp
    hip/TestHIP_UniqueToken.cpp
    hip/TestHIP_ViewAPI_a.cpp
    hip/TestHIP_ViewAPI_b.cpp
    hip/TestHIP_ViewAPI_e.cpp
    hip/TestHIP_ViewLayoutStrideAssignment.cpp
    hip/TestHIP_WorkGraph.cpp
  )

  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HIP
    SOURCES
      ${HIP_SOURCES}
      hip/TestHIPHostPinned_ViewAPI_a.cpp
      hip/TestHIPHostPinned_ViewAPI_b.cpp
      hip/TestHIPHostPinned_ViewAPI_c.cpp
      hip/TestHIPHostPinned_ViewAPI_d.cpp
      hip/TestHIPHostPinned_ViewAPI_e.cpp
      hip/TestHIPHostPinned_ViewCopy.cpp
      hip/TestHIPHostPinned_ViewMapping_a.cpp
      hip/TestHIPHostPinned_ViewMapping_b.cpp
      hip/TestHIPHostPinned_ViewMapping_subview.cpp
  )
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_HIPInterOpInit
    SOURCES
      UnitTestMain.cpp
      hip/TestHIP_InterOp_Init.cpp
  )
endif()

SET(DEFAULT_DEVICE_SOURCES
  UnitTestMainInit.cpp
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
)

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Default
  SOURCES ${DEFAULT_DEVICE_SOURCES}
)

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_PushFinalizeHook
  SOURCES
    UnitTest_PushFinalizeHook.cpp
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

foreach(INITTESTS_NUM RANGE 1 16)
KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DefaultInit_${INITTESTS_NUM}
  SOURCES UnitTestMain.cpp default/TestDefaultDeviceTypeInit_${INITTESTS_NUM}.cpp
)
endforeach(INITTESTS_NUM)

if (KOKKOS_ENABLE_HWLOC)
KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_HWLOC
  SOURCES UnitTestMain.cpp  TestHWLOC.cpp
)
endif()

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_HostBarrier
  SOURCES UnitTestMain.cpp  TestHostBarrier.cpp
)

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

  TARGET_INCLUDE_DIRECTORIES( ${PACKAGE_NAME}_IncrementalTest_${DEVICE} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/incremental )

ENDFUNCTION()

FOREACH (DEVICE ${KOKKOS_ENABLED_DEVICES})
  KOKKOS_ADD_INCREMENTAL_TEST(${DEVICE})
ENDFOREACH()

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_CTestDevice
  SOURCES UnitTestMain.cpp  TestCTestDevice.cpp
)
