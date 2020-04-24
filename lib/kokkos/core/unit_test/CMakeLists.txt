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
KOKKOS_TARGET_COMPILE_DEFINITIONS(kokkos_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0")
KOKKOS_TARGET_INCLUDE_DIRECTORIES(kokkos_gtest PUBLIC ${GTEST_SOURCE_DIR})


#
# Define the tests
#

#I will leave these alone for now because I don't need transitive dependencies on tests
KOKKOS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

foreach(Tag Threads;Serial;OpenMP;Cuda;HPX)
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

if(Kokkos_ENABLE_QTHREADS)
  KOKKOS_ADD_EXECUTABLE_AND_TEST(
    UnitTest_Qthreads
    SOURCES
      UnitTestMainInit.cpp
      qthreads/TestQthreads_Atomics.cpp
      qthreads/TestQthreads_Complex.cpp
      qthreads/TestQthreads_DeepCopyAlignment.cpp
      qthreads/TestQthreads_Other.cpp
      qthreads/TestQthreads_Reductions.cpp
      qthreads/TestQthreads_Reducers_a.cpp
      qthreads/TestQthreads_Reducers_b.cpp
      qthreads/TestQthreads_Reducers_c.cpp
      qthreads/TestQthreads_Reducers_d.cpp
      qthreads/TestQthreads_SubView_a.cpp
      qthreads/TestQthreads_SubView_b.cpp
      qthreads/TestQthreads_SubView_c01.cpp
      qthreads/TestQthreads_SubView_c02.cpp
      qthreads/TestQthreads_SubView_c03.cpp
      qthreads/TestQthreads_SubView_c04.cpp
      qthreads/TestQthreads_SubView_c05.cpp
      qthreads/TestQthreads_SubView_c06.cpp
      qthreads/TestQthreads_SubView_c07.cpp
      qthreads/TestQthreads_SubView_c08.cpp
      qthreads/TestQthreads_SubView_c09.cpp
      qthreads/TestQthreads_SubView_c10.cpp
      qthreads/TestQthreads_SubView_c11.cpp
      qthreads/TestQthreads_SubView_c12.cpp
      qthreads/TestQthreads_SubView_c13.cpp
      qthreads/TestQthreads_Team.cpp
      qthreads/TestQthreads_View_64bit.cpp
      qthreads/TestQthreads_ViewAPI_a.cpp
      qthreads/TestQthreads_ViewAPI_b.cpp
      qthreads/TestQthreads_ViewAPI_c.cpp
      qthreads/TestQthreads_ViewAPI_d.cpp
      qthreads/TestQthreads_ViewAPI_e.cpp
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

KOKKOS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Default
  SOURCES
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
if(NOT KOKKOS_HAS_TRILINOS)
SET_PROPERTY(TARGET StackTraceTestExec PROPERTY ENABLE_EXPORTS 1)

KOKKOS_ADD_TEST( NAME UnitTest_StackTraceTest_normal
                 EXE  StackTraceTestExec 
                 FAIL_REGULAR_EXPRESSION "FAILED"
               )

KOKKOS_ADD_TEST( NAME UnitTest_StackTraceTest_terminate
                 EXE  StackTraceTestExec
                 FAIL_REGULAR_EXPRESSION "FAILED"
               )

KOKKOS_ADD_TEST( NAME UnitTest_StackTraceTest_generic_term
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



