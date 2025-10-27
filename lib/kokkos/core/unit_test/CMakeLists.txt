#
# Add test-only library for gtest to be reused by all the subpackages
#

if(NOT GTest_FOUND) # fallback to internal gtest
  set(GTEST_SOURCE_DIR ${Kokkos_SOURCE_DIR}/tpls/gtest)

  kokkos_add_test_library(
    kokkos_gtest HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
  )

  target_include_directories(kokkos_gtest SYSTEM PUBLIC ${GTEST_SOURCE_DIR})
  if((NOT (Kokkos_ENABLE_CUDA AND WIN32)) AND (NOT ("${KOKKOS_CXX_COMPILER_ID}" STREQUAL "Fujitsu")))
    target_compile_features(kokkos_gtest PUBLIC cxx_std_14)
  endif()

  # Suppress clang-tidy diagnostics on code that we do not have control over
  if(CMAKE_CXX_CLANG_TIDY)
    set_target_properties(kokkos_gtest PROPERTIES CXX_CLANG_TIDY "")
  endif()

  find_package(Threads QUIET)
  if(TARGET Threads::Threads)
    set_target_properties(kokkos_gtest PROPERTIES INTERFACE_LINK_LIBRARIES Threads::Threads)
  endif()
endif()

#
# Define Incremental Testing Feature Levels
# Define Device name mappings (i.e. what comes after Kokkos:: for the ExecSpace)
#

set(KOKKOS_CUDA_FEATURE_LEVEL 999)
set(KOKKOS_CUDA_NAME Cuda)
set(KOKKOS_HIP_FEATURE_LEVEL 999)
set(KOKKOS_HIP_NAME HIP)
set(KOKKOS_HPX_FEATURE_LEVEL 999)
set(KOKKOS_HPX_NAME Experimental::HPX)
set(KOKKOS_OPENMP_FEATURE_LEVEL 999)
set(KOKKOS_OPENMP_NAME OpenMP)

# FIXME_OPENMPTARGET - The NVIDIA HPC compiler nvc++ only compiles the first 8 incremental tests for the OpenMPTarget backend.
# FIXME_OPENMPTARGET - Clang version 17 fails to compile incremental tests past 12 with verion 17. There is PR for this in upstream already. So it should be fixed by version 18.
if(KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  set(KOKKOS_OPENMPTARGET_FEATURE_LEVEL 10)
elseif(KOKKOS_CXX_COMPILER_ID STREQUAL "Clang" AND KOKKOS_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 17.0.0)
  set(KOKKOS_OPENMPTARGET_FEATURE_LEVEL 12)
else()
  set(KOKKOS_OPENMPTARGET_FEATURE_LEVEL 14)
endif()

set(KOKKOS_OPENMPTARGET_NAME Experimental::OpenMPTarget)
set(KOKKOS_SERIAL_FEATURE_LEVEL 999)
set(KOKKOS_SERIAL_NAME Serial)
set(KOKKOS_SYCL_FEATURE_LEVEL 999)
set(KOKKOS_SYCL_NAME SYCL)
set(KOKKOS_THREADS_FEATURE_LEVEL 999)
set(KOKKOS_THREADS_NAME Threads)
# FIXME_OPENACC - The Clang compiler only compiles the first 9 incremental tests for the OpenACC backend.
if(KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
  set(KOKKOS_OPENACC_FEATURE_LEVEL 9)
else()
  set(KOKKOS_OPENACC_FEATURE_LEVEL 17)
endif()

set(KOKKOS_OPENACC_NAME Experimental::OpenACC)

#
# Define the tests
#

#I will leave these alone for now because I don't need transitive dependencies on tests
kokkos_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkos_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
kokkos_include_directories(${KOKKOS_SOURCE_DIR}/core/unit_test/category_files)

set(COMPILE_ONLY_SOURCES
    TestArray.cpp
    TestCreateMirror.cpp
    TestDetectionIdiom.cpp
    TestBitManipulation.cpp
    TestInterOp.cpp
    TestRangePolicyCTAD.cpp
    TestStringManipulation.cpp
    TestVersionMacros.cpp
    TestViewRank.cpp
    TestViewTypeTraits.cpp
    TestViewTypedefs.cpp
    TestTypeInfo.cpp
    TestTypeList.cpp
    TestMDRangePolicyCTAD.cpp
    TestTeamPolicyCTAD.cpp
    TestTeamMDRangePolicyCTAD.cpp
    TestNestedReducerCTAD.cpp
    view/TestBasicViewMDSpanConversion.cpp
    view/TestViewCustomization.cpp
    view/TestExtentsDatatypeConversion.cpp
)

if(NOT Kokkos_ENABLE_IMPL_MDSPAN OR KOKKOS_CXX_COMPILER_ID STREQUAL "Intel")
  list(REMOVE_ITEM COMPILE_ONLY_SOURCES view/TestBasicViewMDSpanConversion.cpp)
endif()

if(NOT Kokkos_ENABLE_IMPL_MDSPAN OR Kokkos_ENABLE_IMPL_VIEW_LEGACY)
  list(REMOVE_ITEM COMPILE_ONLY_SOURCES view/TestViewCustomization.cpp)
endif()

#testing if windows.h and Kokkos_Core.hpp can be included
if(WIN32)
  list(APPEND COMPILE_ONLY_SOURCES TestWindowsInclude.cpp)
endif()

if(Kokkos_INSTALL_TESTING) # FIXME Kokkos_ and KOKKOS_ variables are out of sync
  if((Kokkos_CXX_COMPILER_ID STREQUAL "Intel" AND Kokkos_CXX_COMPILER_VERSION VERSION_LESS 2021.2.0)
     OR (Kokkos_CXX_COMPILER_ID STREQUAL "NVIDIA" AND Kokkos_CXX_COMPILER_VERSION VERSION_LESS 11.3.0)
     OR (Kokkos_CXX_COMPILER_ID STREQUAL "NVIDIA" AND Kokkos_CXX_HOST_COMPILER_ID STREQUAL "MSVC")
  )
    list(REMOVE_ITEM COMPILE_ONLY_SOURCES TestTypeInfo.cpp)
  endif()
else()
  if((KOKKOS_CXX_COMPILER_ID STREQUAL "Intel" AND KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 2021.2.0)
     OR (KOKKOS_CXX_COMPILER_ID STREQUAL "NVIDIA" AND KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 11.3.0)
     OR (KOKKOS_CXX_COMPILER_ID STREQUAL "NVIDIA" AND KOKKOS_CXX_HOST_COMPILER_ID STREQUAL "MSVC")
  )
    list(REMOVE_ITEM COMPILE_ONLY_SOURCES TestTypeInfo.cpp)
  endif()
endif()

if(Kokkos_ENABLE_OPENMPTARGET)
  list(REMOVE_ITEM COMPILE_ONLY_SOURCES TestNestedReducerCTAD.cpp)
endif()
kokkos_add_executable(CoreTestCompileOnly SOURCES TestCompileMain.cpp ${COMPILE_ONLY_SOURCES})

foreach(Tag Threads;Serial;OpenMP;Cuda;HPX;OpenMPTarget;OpenACC;HIP;SYCL)
  string(TOUPPER ${Tag} DEVICE)
  string(TOLOWER ${Tag} dir)

  if(Kokkos_ENABLE_${DEVICE})
    set(dir ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    file(MAKE_DIRECTORY ${dir})
    # Needed to split this for Windows NVCC, since it ends up putting everything on the
    # command line in an intermediate compilation step even if CMake generated a response
    # file. That then exceeded the shell command line max length.
    set(${Tag}_SOURCES1A)
    foreach(
      Name
      Abort
      ArrayOps
      AtomicOperations_complexdouble
      AtomicOperations_complexfloat
      AtomicOperations_double
      AtomicOperations_float
      AtomicOperations_int
      AtomicOperations_longint
      AtomicOperations_longlongint
      AtomicOperations_shared
      AtomicOperations_unsignedint
      AtomicOperations_unsignedlongint
      Atomics
      AtomicViews
      BitManipulationBuiltins
      BlockSizeDeduction
      CheckedIntegerOps
      CommonPolicyConstructors
      CommonPolicyInterface
      Complex
      Concepts
      Crs
      DeepCopyAlignment
      ExecSpacePartitioning
      ExecSpaceThreadSafety
      ExecutionSpace
      FunctorAnalysis
      Graph
      HostSharedPtr
      HostSharedPtrAccessOnDevice
      JoinBackwardCompatibility
      LocalDeepCopy
      MathematicalConstants
      MathematicalFunctions1
      MathematicalFunctions2
      MathematicalFunctions3
      MathematicalSpecialFunctions
    )
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp "#include <Test${Tag}_Category.hpp>\n" "#include <Test${Name}.hpp>\n")
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES1A ${file})
    endforeach()

    set(${Tag}_SOURCES1B)
    set(${Tag}_TESTNAMES1B
        MDRange_a
        MDRange_b
        MDRange_c
        MDRange_d
        MDRange_e
        MDRange_f
        MDRange_g
        MDRangePolicyConstructors
        MDRangeReduce
        MDSpan
        MDSpanAtomicAccessor
        MDSpanConversion
        MinMaxClamp
        NumericTraits
        OccupancyControlTrait
        Other
        ParallelScanRangePolicy
        Printf
        QuadPrecisionMath
        RangePolicy
        RangePolicyConstructors
        RangePolicyRequire
        ReducerCTADs
        Reducers_a
        Reducers_b
        Reducers_c
        Reducers_d
        Reducers_e
        Reductions
        Reductions_DeviceView
        SharedAlloc
        SpaceAwareAccessorAccessViolation
        SpaceAwareAccessor
        Swap
    )
    if(NOT Kokkos_ENABLE_IMPL_MDSPAN)
      list(REMOVE_ITEM ${Tag}_TESTNAMES1B MDSpanAtomicAccessor MDSpanConversion SpaceAwareAccessorAccessViolation
           SpaceAwareAccessor
      )
    endif()
    # This test case causes MSVC to fail with "number of sections exceeded object file format limit"
    if(MSVC)
      list(REMOVE_ITEM ${Tag}_TESTNAMES1B Reducers_d)
    endif()
    foreach(Name IN LISTS ${Tag}_TESTNAMES1B)
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp "#include <Test${Tag}_Category.hpp>\n" "#include <Test${Name}.hpp>\n")
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES1B ${file})
    endforeach()

    set(${Tag}_SOURCES2A)
    set(${Tag}_TESTNAMES2A
        TeamBasic
        TeamCombinedReducers
        TeamMDRange
        TeamPolicyConstructors
        TeamReductionScan
        TeamScan
        TeamScratch
        TeamTeamSize
        TeamVectorRange
        UniqueToken
        View_64bit
        ViewAPI_a
        ViewAPI_b
        ViewAPI_c
        ViewAPI_d
        ViewAPI_e
        ViewBadAlloc
        ViewCopy_a
        ViewCopy_b
        ViewCopy_c
        ViewCtorDimMatch
        ViewCtorProp
        ViewEmptyRuntimeUnmanaged
        ViewHooks
        ViewLayoutStrideAssignment
        ViewMapping_a
        ViewMapping_b
        ViewMapping_subview
        ViewMemoryAccessViolation
        ViewOfClass
        ViewOfViews
        ViewOutOfBoundsAccess
        ViewResize
        WorkGraph
        WithoutInitializing
    )
    # Workaround to internal compiler error with intel classic compilers
    # when using -no-ip flag in ViewCopy_c
    # See issue: https://github.com/kokkos/kokkos/issues/7084
    if(KOKKOS_CXX_COMPILER_ID STREQUAL Intel)
      list(REMOVE_ITEM ${Tag}_TESTNAMES2A ViewCopy_c)
    endif()
    # FIXME: we will eventually likely remove ViewHooks entirely, but for now
    # ViewHooks are not compatible
    if(NOT Kokkos_ENABLE_IMPL_VIEW_LEGACY)
      list(REMOVE_ITEM ${Tag}_TESTNAMES2A ViewHooks)
    endif()
    foreach(Name IN LISTS ${Tag}_TESTNAMES2A)
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp "#include <Test${Tag}_Category.hpp>\n" "#include <Test${Name}.hpp>\n")
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES2A ${file})
    endforeach()

    set(TagHostAccessible ${Tag})
    if(Tag STREQUAL "Cuda")
      set(TagHostAccessible CudaUVM)
    elseif(Tag STREQUAL "HIP")
      set(TagHostAccessible HIPManaged)
    elseif(Tag STREQUAL "SYCL")
      set(TagHostAccessible SYCLSharedUSM)
    endif()

    set(${Tag}_SOURCES2B)
    foreach(
      Name
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
      file(WRITE ${dir}/dummy.cpp "#include <Test${TagHostAccessible}_Category.hpp>\n" "#include <Test${Name}.hpp>\n")
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES2B ${file})
    endforeach()

    set(${Tag}_SOURCES2C)
    foreach(Name SubView_c06 SubView_c07 SubView_c08 SubView_c09)
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp "#include <Test${TagHostAccessible}_Category.hpp>\n" "#include <Test${Name}.hpp>\n")
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES2C ${file})
    endforeach()

    set(${Tag}_SOURCES2D)
    foreach(Name SubView_c10 SubView_c11 SubView_c12 SubView_c13 SubView_c14 SubView_c15)
      set(file ${dir}/Test${Tag}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp "#include <Test${TagHostAccessible}_Category.hpp>\n" "#include <Test${Name}.hpp>\n")
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${Tag}_SOURCES2D ${file})
    endforeach()

    set(${Tag}_SOURCES1 ${${Tag}_SOURCES1A} ${${Tag}_SOURCES1B})
    set(${Tag}_SOURCES2 ${${Tag}_SOURCES2A} ${${Tag}_SOURCES2B} ${${Tag}_SOURCES2C} ${${Tag}_SOURCES2D})
    set(${Tag}_SOURCES ${${Tag}_SOURCES1} ${${Tag}_SOURCES2})

    # ViewSupport should eventually contain the new implementation
    # detail tests for the mdspan based View
    if(Kokkos_ENABLE_IMPL_MDSPAN)
      set(BV_TestNames BasicView ReferenceCountedAccessor ReferenceCountedDataHandle)
      if(NOT Kokkos_ENABLE_IMPL_VIEW_LEGACY AND NOT Kokkos_ENABLE_OPENACC)
        list(APPEND BV_TestNames ViewCustomizationAccessorArg ViewCustomizationAllocationType
             ViewCustomizationAccessorFromMapping
        )
      endif()
      set(${Tag}_VIEWSUPPORT)
      foreach(Name IN LISTS BV_TestNames)
        set(file ${dir}/Test${Tag}_View_${Name}.cpp)
        # Write to a temporary intermediate file and call configure_file to avoid
        # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
        file(WRITE ${dir}/dummy.cpp "#include <Test${Tag}_Category.hpp>\n" "#include <view/Test${Name}.hpp>\n")
        configure_file(${dir}/dummy.cpp ${file})
        list(APPEND ${Tag}_VIEWSUPPORT ${file})
      endforeach()
      kokkos_add_executable_and_test(CoreUnitTest_${Tag}_ViewSupport SOURCES UnitTestMainInit.cpp ${${Tag}_VIEWSUPPORT})
    endif()
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
    foreach(
      Name
      SharedAlloc
      ViewAPI_a
      ViewAPI_b
      ViewAPI_c
      ViewAPI_d
      ViewAPI_e
      ViewCopy_a
      ViewCopy_b
      ViewCopy_c
      ViewMapping_a
      ViewMapping_b
      ViewMapping_subview
    )
      set(file ${dir}/Test${DEVICE}${SPACE}_${Name}.cpp)
      # Write to a temporary intermediate file and call configure_file to avoid
      # updating timestamps triggering unnecessary rebuilds on subsequent cmake runs.
      file(WRITE ${dir}/dummy.cpp "#include <Test${DEVICE}${SPACE}_Category.hpp>\n" "#include <Test${Name}.hpp>\n")
      configure_file(${dir}/dummy.cpp ${file})
      list(APPEND ${DEVICE}_SOURCES3 ${file})
    endforeach()
    list(APPEND ${DEVICE}_SOURCES ${${DEVICE}_SOURCES3})
  endif()
endforeach()

# Disable non-compiling tests based on clang version.
if(Kokkos_ENABLE_OPENMPTARGET)
  list(
    REMOVE_ITEM
    OpenMPTarget_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Other.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamCombinedReducers.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamReductionScan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_WorkGraph.cpp
    IF
    (KOKKOS_CXX_COMPILER_ID STREQUAL "Clang" AND KOKKOS_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 16.0.0)
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_shared.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MinMaxClamp.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_LocalDeepCopy.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_d.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamVectorRange.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_e.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamScan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamBasic.cpp
    IF
    (KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 17.0.3)
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c01.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c02.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SubView_c03.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_d.cpp
    endif
    ()
    endif
    ()
    # FIXME_OPENMPTARGET_CRAY: The following tests fail at compile time when the OpenMPTarget backend is enabled with the Cray compiler.
    # Atomic compare/exchange is used in these tests which can be one of the reasons for the compilation failures.
    IF
    (KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_complexdouble.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_UniqueToken.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_SharedAlloc.cpp
    ENDIF
    ()
  )
endif()

# FIXME_OPENMPTARGET - MinMaxClamp fails even with the host backend when OpenMPTarget backend is enabled.
# FIXME_OPENMPTARGET - Unsure of the reason as of now.
if(KOKKOS_CXX_COMPILER_ID STREQUAL "Clang" AND KOKKOS_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 16.0.0)
  if(Kokkos_ENABLE_OPENMPTARGET AND Kokkos_ENABLE_OPENMP)
    list(REMOVE_ITEM OpenMP_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/openmp/TestOpenMP_MinMaxClamp.cpp)
  endif()
  if(Kokkos_ENABLE_OPENMPTARGET AND Kokkos_ENABLE_SERIAL)
    list(REMOVE_ITEM Serial_SOURCES1 ${CMAKE_CURRENT_BINARY_DIR}/serial/TestSerial_MinMaxClamp.cpp)
  endif()
endif()

if(Kokkos_ENABLE_OPENACC)
  list(
    REMOVE_ITEM
    OpenACC_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_complexdouble.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_complexfloat.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Crs.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_JoinBackwardCompatibility.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_LocalDeepCopy.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Other.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamCombinedReducers.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamMDRange.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamReductionScan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamScan.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamVectorRange.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewAPI_e.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewMapping_subview.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewOfClass.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_WorkGraph.cpp
  )
endif()

# FIXME_OPENMPTARGET - Comment non-passing tests with amdclang++
# FIXME_OPENMPTARGET - Need to check on GFX1030 and GFX1100 architectures
if(KOKKOS_ARCH_VEGA)
  set(KOKKOS_AMDGPU_ARCH TRUE)
endif()
if(KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL Clang AND KOKKOS_AMDGPU_ARCH)
  list(REMOVE_ITEM OpenMPTarget_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_e.cpp)
endif()

# FIXME_OPENMPTARGET - Comment non-passing tests with the NVIDIA HPC compiler nvc++
if(KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  list(
    REMOVE_ITEM
    OpenMPTarget_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/default/TestDefaultDeviceType_a1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/default/TestDefaultDeviceType_b1.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_double.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_float.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_int.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_longint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_longlongint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_unsignedint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicOperations_unsignedlongint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Atomics.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_AtomicViews.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_BlockSizeDeduction.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_CommonPolicyConstructors.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_CommonPolicyInterface.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_DeepCopyAlignment.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_HostSharedPtr.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_HostSharedPtrAccessOnDevice.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MathematicalFunctions.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MDRange_a.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MDRange_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MDRange_c.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MDRange_d.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_MDRangePolicyConstructors.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_NumericTraits.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_RangePolicy.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_RangePolicyConstructors.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_RangePolicyRequire.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_a.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_c.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_d.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reducers_e.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_Reductions_DeviceView.cpp
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
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamBasic.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamPolicyConstructors.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamScratch.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamTeamSize.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_TeamVectorRange.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_UniqueToken.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_a.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_c.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_d.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewAPI_f.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewMapping_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openmptarget/TestOpenMPTarget_ViewResize.cpp
  )
endif()

# FIXME_OPENACC - Comment non-passing tests with the NVIDIA HPC compiler nvc++
if(KOKKOS_ENABLE_OPENACC AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  list(
    REMOVE_ITEM
    OpenACC_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/default/TestDefaultDeviceType_a1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/default/TestDefaultDeviceType_b1.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_shared.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_BlockSizeDeduction.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_DeepCopyAlignment.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_HostSharedPtr.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_HostSharedPtrAccessOnDevice.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MathematicalFunctions1.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MathematicalFunctions2.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MathematicalFunctions3.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MDRange_c.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MDRange_f.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_NumericTraits.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_RangePolicy.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_RangePolicyRequire.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Reducers_a.cpp #fails if NVHPC V22.5 or lower.
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Reducers_d.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Reductions.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Reductions_DeviceView.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamBasic.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamScratch.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamTeamSize.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_UniqueToken.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewResize.cpp
  )
endif()

# FIXME_OPENACC - Comment non-passing tests with the Clang compiler
if(KOKKOS_ENABLE_OPENACC AND KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
  list(
    REMOVE_ITEM
    OpenACC_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/default/TestDefaultDeviceType_a1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/default/TestDefaultDeviceType_b1.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_double.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_float.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_int.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_longint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_longlongint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_shared.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_unsignedint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicOperations_unsignedlongint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Atomics.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_AtomicViews.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_BlockSizeDeduction.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_DeepCopyAlignment.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_HostSharedPtrAccessOnDevice.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MathematicalFunctions1.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MathematicalFunctions2.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MDRange_c.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MDRange_f.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_NumericTraits.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_RangePolicy.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_RangePolicyRequire.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Reducers_a.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Reducers_c.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Reducers_d.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Reductions.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Reductions_DeviceView.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamBasic.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamScratch.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_TeamTeamSize.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_UniqueToken.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewMapping_b.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewResize.cpp
    # This test is not removed above for OpenACC+NVHPC but all its TEST
    # functions are not compiled for the case of KOKKOS_COMPILER_NVHPC.
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewCtorDimMatch.cpp
    # These tests are not removed above for OpenACC+NVHPC.
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Abort.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_Complex.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ExecutionSpace.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ExecSpacePartitioning.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MathematicalConstants.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MathematicalSpecialFunctions.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MinMaxClamp.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewLayoutStrideAssignment.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewMapping_a.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewMemoryAccessViolation.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_WithoutInitializing.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewAPI_d.cpp
    #Below test is disabled because it uses atomic operations not supported by Clacc.
    ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ExecSpaceThreadSafety.cpp
  )
  # When tested on a systme with AMD MI60 GPU and ROCm V5.4.0, these cause
  # clang-linker-wrapper to hang for a long time while building the unit tests.
  # In some cases, including them caused the build not to complete after an hour,
  # but excluding them permitted the build to finish in 1.5 mins or less.
  if(KOKKOS_AMDGPU_ARCH)
    list(
      REMOVE_ITEM
      OpenACC_SOURCES
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_BitManipulationBuiltins.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_MathematicalFunctions3.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ParallelScanRangePolicy.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_SubView_c04.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_SubView_c05.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_SubView_c06.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_SubView_c07.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_SubView_c08.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_SubView_c09.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_SubView_c10.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_SubView_c11.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_SubView_c12.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewAPI_b.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/openacc/TestOpenACC_ViewAPI_c.cpp
    )
  endif()
  # Fails serial.atomics_tpetra_max_abs when we test with Clacc.
  list(REMOVE_ITEM Serial_SOURCES1 ${CMAKE_CURRENT_BINARY_DIR}/serial/TestSerial_Atomics.cpp)
endif()

if(Kokkos_ENABLE_SERIAL)
  if(Kokkos_ENABLE_DEPRECATED_CODE_4)
    set(Serial_EXTRA_SOURCES serial/TestSerial_Task.cpp)
  else()
    set(Serial_EXTRA_SOURCES)
  endif()

  kokkos_add_executable_and_test(
    CoreUnitTest_Serial1 SOURCES UnitTestMainInit.cpp ${Serial_SOURCES1} ${Serial_EXTRA_SOURCES}
  )
  kokkos_add_executable_and_test(CoreUnitTest_Serial2 SOURCES UnitTestMainInit.cpp ${Serial_SOURCES2})
endif()

if(Kokkos_ENABLE_THREADS)
  kokkos_add_executable_and_test(CoreUnitTest_Threads SOURCES ${Threads_SOURCES} UnitTestMainInit.cpp)
endif()

if(Kokkos_ENABLE_OPENMP)
  if(Kokkos_ENABLE_DEPRECATED_CODE_4)
    set(OpenMP_EXTRA_SOURCES openmp/TestOpenMP_Task.cpp)
  else()
    set(OpenMP_EXTRA_SOURCES)
  endif()
  kokkos_add_executable_and_test(
    CoreUnitTest_OpenMP SOURCES UnitTestMainInit.cpp ${OpenMP_SOURCES} ${OpenMP_EXTRA_SOURCES}
  )
  kokkos_add_executable_and_test(CoreUnitTest_OpenMPInterOp SOURCES UnitTestMain.cpp openmp/TestOpenMP_InterOp.cpp)
endif()

if(Kokkos_ENABLE_HPX)
  if(Kokkos_ENABLE_DEPRECATED_CODE_4)
    set(HPX_EXTRA_SOURCES hpx/TestHPX_Task.cpp)
  else()
    set(HPX_EXTRA_SOURCES)
  endif()

  kokkos_add_executable_and_test(CoreUnitTest_HPX SOURCES UnitTestMainInit.cpp ${HPX_SOURCES} ${HPX_EXTRA_SOURCES})
  kokkos_add_executable_and_test(CoreUnitTest_HPXInterOp SOURCES UnitTestMain.cpp hpx/TestHPX_InterOp.cpp)
  kokkos_add_executable_and_test(
    CoreUnitTest_HPX_IndependentInstances
    SOURCES
    UnitTestMainInit.cpp
    hpx/TestHPX_IndependentInstances.cpp
    hpx/TestHPX_IndependentInstancesDelayedExecution.cpp
    hpx/TestHPX_IndependentInstancesInstanceIds.cpp
    hpx/TestHPX_IndependentInstancesRefCounting.cpp
    hpx/TestHPX_IndependentInstancesSynchronization.cpp
  )
  if(Kokkos_ENABLE_DEPRECATED_CODE_4)
    kokkos_add_executable_and_test(CoreUnitTest_HPX_InParallel SOURCES UnitTestMainInit.cpp hpx/TestHPX_InParallel.cpp)
  endif()
endif()

if(Kokkos_ENABLE_OPENMPTARGET)
  kokkos_add_executable_and_test(CoreUnitTest_OpenMPTarget SOURCES UnitTestMainInit.cpp ${OpenMPTarget_SOURCES})
endif()

if(Kokkos_ENABLE_OPENACC)
  kokkos_add_executable_and_test(CoreUnitTest_OpenACC SOURCES UnitTestMainInit.cpp ${OpenACC_SOURCES})
endif()

if(Kokkos_ENABLE_CUDA)
  kokkos_add_executable_and_test(
    CoreUnitTest_Cuda1 SOURCES UnitTestMainInit.cpp ${Cuda_SOURCES1} cuda/TestCuda_ReducerViewSizeLimit.cpp
  )

  kokkos_add_executable_and_test(CoreUnitTest_Cuda2 SOURCES UnitTestMainInit.cpp ${Cuda_SOURCES2})

  if(Kokkos_ENABLE_DEPRECATED_CODE_4)
    set(Cuda_EXTRA_SOURCES cuda/TestCuda_Task.cpp)
  else()
    set(Cuda_EXTRA_SOURCES)
  endif()

  kokkos_add_executable_and_test(
    CoreUnitTest_Cuda3
    SOURCES
    UnitTestMainInit.cpp
    ${Cuda_EXTRA_SOURCES}
    cuda/TestCuda_TeamScratchStreams.cpp
    ${Cuda_SOURCES3}
    cuda/TestCuda_Spaces.cpp
    ${Cuda_SOURCES_SHAREDSPACE}
  )

  # This test seeks to make sure that `desul::ensure_cuda_lock_arrays_on_device` is called before a graph submission.
  # It is in a separate file to ensure that the test takes place in a separate compilation unit, in which this function
  # has not been called yet by Kokkos::initialize or by a preceding kernel launch.
  set(file ${CMAKE_CURRENT_BINARY_DIR}/cuda/TestCuda_GraphAtomicLocks.cpp)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cuda/dummy.cpp "#include <TestCuda_Category.hpp>\n"
                                                        "#include <TestGraphAtomicLocks.hpp>\n"
  )
  configure_file(${CMAKE_CURRENT_BINARY_DIR}/cuda/dummy.cpp ${file})
  kokkos_add_executable_and_test(CoreUnitTest_CudaGraphAtomicLocks SOURCES UnitTestMainInit.cpp ${file})

  kokkos_add_executable_and_test(
    CoreUnitTest_CudaTimingBased SOURCES UnitTestMainInit.cpp cuda/TestCuda_DebugSerialExecution.cpp
    cuda/TestCuda_DebugPinUVMSpace.cpp
  )

  kokkos_add_executable_and_test(CoreUnitTest_CudaInterOpInit SOURCES UnitTestMain.cpp cuda/TestCuda_InterOp_Init.cpp)
  kokkos_add_executable_and_test(
    CoreUnitTest_CudaInterOpStreams SOURCES UnitTestMain.cpp cuda/TestCuda_InterOp_Streams.cpp
  )
  kokkos_add_executable_and_test(
    CoreUnitTest_CudaInterOpStreamsMultiGPU SOURCES UnitTestMainInit.cpp cuda/TestCuda_InterOp_StreamsMultiGPU.cpp
  )

  kokkos_add_executable_and_test(
    CoreUnitTest_CudaInterOpGraph SOURCES UnitTestMainInit.cpp cuda/TestCuda_InterOp_Graph.cpp
  )
endif()

if(Kokkos_ENABLE_HIP)
  kokkos_add_executable_and_test(
    CoreUnitTest_HIP
    SOURCES
    UnitTestMainInit.cpp
    ${HIP_SOURCES}
    hip/TestHIP_ScanUnit.cpp
    hip/TestHIP_SharedResourceLock.cpp
    hip/TestHIP_Spaces.cpp
    hip/TestHIP_Memory_Requirements.cpp
    hip/TestHIP_TeamScratchStreams.cpp
    hip/TestHIP_AsyncLauncher.cpp
    hip/TestHIP_BlocksizeDeduction.cpp
    hip/TestHIP_UnifiedMemory_ZeroMemset.cpp
  )

  set(file ${CMAKE_CURRENT_BINARY_DIR}/hip/TestHIP_GraphAtomicLocks.cpp)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hip/dummy.cpp "#include <TestHIP_Category.hpp>\n"
                                                       "#include <TestGraphAtomicLocks.hpp>\n"
  )
  configure_file(${CMAKE_CURRENT_BINARY_DIR}/hip/dummy.cpp ${file})
  kokkos_add_executable_and_test(CoreUnitTest_HIPGraphAtomicLocks SOURCES UnitTestMainInit.cpp ${file})

  kokkos_add_executable_and_test(CoreUnitTest_HIPInterOpInit SOURCES UnitTestMain.cpp hip/TestHIP_InterOp_Init.cpp)
  kokkos_add_executable_and_test(
    CoreUnitTest_HIPInterOpStreams SOURCES UnitTestMain.cpp hip/TestHIP_InterOp_Streams.cpp
  )
  kokkos_add_executable_and_test(
    CoreUnitTest_HIPInterOpGraph SOURCES UnitTestMainInit.cpp hip/TestHIP_InterOp_Graph.cpp
  )
  if(KOKKOS_CXX_COMPILER_ID STREQUAL HIPCC AND KOKKOS_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 5.6)
    # MultiGPU only works for rocm 5.6+
    kokkos_add_executable_and_test(
      CoreUnitTest_HIPInterOpStreamsMultiGPU SOURCES UnitTestMainInit.cpp hip/TestHIP_InterOp_StreamsMultiGPU.cpp
    )
  endif()
endif()

if(Kokkos_ENABLE_SYCL)
  list(REMOVE_ITEM SYCL_SOURCES2A ${CMAKE_CURRENT_BINARY_DIR}/sycl/TestSYCL_WorkGraph.cpp)

  kokkos_add_executable_and_test(CoreUnitTest_SYCL1A SOURCES UnitTestMainInit.cpp ${SYCL_SOURCES1A})

  kokkos_add_executable_and_test(CoreUnitTest_SYCL1B SOURCES UnitTestMainInit.cpp ${SYCL_SOURCES1B})

  kokkos_add_executable_and_test(CoreUnitTest_SYCL2A SOURCES UnitTestMainInit.cpp ${SYCL_SOURCES2A})

  kokkos_add_executable_and_test(CoreUnitTest_SYCL2B SOURCES UnitTestMainInit.cpp ${SYCL_SOURCES2B})

  kokkos_add_executable_and_test(CoreUnitTest_SYCL2C SOURCES UnitTestMainInit.cpp ${SYCL_SOURCES2C})

  kokkos_add_executable_and_test(CoreUnitTest_SYCL2D SOURCES UnitTestMainInit.cpp ${SYCL_SOURCES2D})

  kokkos_add_executable_and_test(
    CoreUnitTest_SYCL3 SOURCES UnitTestMainInit.cpp sycl/TestSYCL_TeamScratchStreams.cpp ${SYCL_SOURCES3}
    sycl/TestSYCL_Spaces.cpp
  )

  kokkos_add_executable_and_test(CoreUnitTest_SYCLInterOpInit SOURCES UnitTestMain.cpp sycl/TestSYCL_InterOp_Init.cpp)
  kokkos_add_executable_and_test(
    CoreUnitTest_SYCLInterOpInit_Context SOURCES UnitTestMainInit.cpp sycl/TestSYCL_InterOp_Init_Context.cpp
  )
  kokkos_add_executable_and_test(
    CoreUnitTest_SYCLInterOpStreams SOURCES UnitTestMain.cpp sycl/TestSYCL_InterOp_Streams.cpp
  )
  kokkos_add_executable_and_test(
    CoreUnitTest_SYCLInterOpStreamsMultiGPU SOURCES UnitTestMainInit.cpp sycl/TestSYCL_InterOp_StreamsMultiGPU.cpp
  )

  if(KOKKOS_IMPL_HAVE_SYCL_EXT_ONEAPI_GRAPH AND NOT Kokkos_ENABLE_IMPL_SYCL_OUT_OF_ORDER_QUEUES)
    kokkos_add_executable_and_test(
      CoreUnitTest_SYCLInterOpGraph SOURCES UnitTestMainInit.cpp sycl/TestSYCL_InterOp_Graph.cpp
    )
  endif()
endif()

set(DEFAULT_DEVICE_SOURCES
    UnitTestMainInit.cpp
    TestCStyleMemoryManagement.cpp
    TestSharedSpace.cpp
    TestSharedHostPinnedSpace.cpp
    TestCompilerMacros.cpp
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
    default/TestDefaultDeviceTypeResize.cpp
    default/TestDefaultDeviceTypeViewAPI.cpp
)
# FIXME_OPENMPTARGET and FIXME_OPENACC do not provide a MemorySpace that can be accessed from all ExecSpaces
# FIXME_SYCL clock_tic does not give the correct timings for cloc_tic
if(KOKKOS_ENABLE_OPENACC OR KOKKOS_ENABLE_OPENMPTARGET OR KOKKOS_ENABLE_SYCL)
  list(REMOVE_ITEM DEFAULT_DEVICE_SOURCES TestSharedSpace.cpp)
endif()
# FIXME_OPENMPTARGET and FIXME_OPENACC do not provide a HostPinnedMemorySpace that can be accessed from all ExecSpaces
if(KOKKOS_ENABLE_OPENACC OR KOKKOS_ENABLE_OPENMPTARGET)
  list(REMOVE_ITEM DEFAULT_DEVICE_SOURCES TestSharedHostPinnedSpace.cpp)
endif()

# FIXME_OPENMPTARGET, FIXME_OPENACC - Comment non-passing tests with the NVIDIA HPC compiler nvc++
if((KOKKOS_ENABLE_OPENMPTARGET OR KOKKOS_ENABLE_OPENACC) AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
  list(
    REMOVE_ITEM
    DEFAULT_DEVICE_SOURCES
    default/TestDefaultDeviceType_a1.cpp
    default/TestDefaultDeviceType_b1.cpp
    default/TestDefaultDeviceType_c1.cpp
    default/TestDefaultDeviceType_a2.cpp
    default/TestDefaultDeviceType_b2.cpp
    default/TestDefaultDeviceType_c2.cpp
    default/TestDefaultDeviceType_a3.cpp
    default/TestDefaultDeviceType_b3.cpp
    default/TestDefaultDeviceType_c3.cpp
    default/TestDefaultDeviceTypeResize.cpp
    default/TestDefaultDeviceTypeViewAPI.cpp
  )
endif()

# FIXME_OPENACC - Comment non-passing tests with the Clang compiler
if(KOKKOS_ENABLE_OPENACC AND KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
  list(
    REMOVE_ITEM
    DEFAULT_DEVICE_SOURCES
    default/TestDefaultDeviceType_a1.cpp
    default/TestDefaultDeviceType_b1.cpp
    default/TestDefaultDeviceType_c1.cpp
    default/TestDefaultDeviceType_a2.cpp
    default/TestDefaultDeviceType_b2.cpp
    default/TestDefaultDeviceType_c2.cpp
    default/TestDefaultDeviceType_a3.cpp
    default/TestDefaultDeviceType_b3.cpp
    default/TestDefaultDeviceType_c3.cpp
    default/TestDefaultDeviceTypeResize.cpp
    default/TestDefaultDeviceTypeViewAPI.cpp
  )
endif()

kokkos_add_executable_and_test(CoreUnitTest_Default SOURCES ${DEFAULT_DEVICE_SOURCES})

kokkos_add_executable_and_test(
  CoreUnitTest_InitializeFinalize
  SOURCES
  UnitTestMain.cpp
  TestExecutionEnvironmentNonInitializedOrFinalized.cpp
  TestInitializationSettings.cpp
  TestInitializeFinalize.cpp
  TestKokkosHelpCausesNormalProgramTermination.cpp
  TestLegionInitialization.cpp
  TestParseCmdLineArgsAndEnvVars.cpp
  TestPushFinalizeHook.cpp
  TestScopeGuard.cpp
)

# This test is intended for development and debugging by putting code
# into TestDefaultDeviceDevelop.cpp. By default its empty.
kokkos_add_executable_and_test(CoreUnitTest_Develop SOURCES UnitTestMainInit.cpp default/TestDefaultDeviceDevelop.cpp)

if(KOKKOS_ENABLE_TUNING)
  kokkos_add_executable_and_test(CoreUnitTest_TuningBuiltins SOURCES tools/TestBuiltinTuners.cpp)
  kokkos_add_executable_and_test(CoreUnitTest_TuningBasics SOURCES tools/TestTuning.cpp)
  kokkos_add_executable_and_test(CoreUnitTest_CategoricalTuner SOURCES tools/TestCategoricalTuner.cpp)
endif()

set(KOKKOSP_SOURCES UnitTestMainInit.cpp tools/TestEventCorrectness.cpp tools/TestKernelNames.cpp
                    tools/TestProfilingSection.cpp tools/TestScopedRegion.cpp tools/TestWithoutInitializing.cpp
)

kokkos_add_executable_and_test(CoreUnitTest_KokkosP SOURCES ${KOKKOSP_SOURCES})
if(KOKKOS_ENABLE_LIBDL)
  kokkos_add_executable_and_test(CoreUnitTest_ToolIndependence SOURCES tools/TestIndependence.cpp)
  target_compile_definitions(Kokkos_CoreUnitTest_ToolIndependence PUBLIC KOKKOS_TOOLS_INDEPENDENT_BUILD)
  kokkos_add_test_library(kokkosprinter-tool SHARED SOURCES tools/printing-tool.cpp)

  if((NOT (Kokkos_ENABLE_CUDA AND WIN32)) AND (NOT ("${KOKKOS_CXX_COMPILER_ID}" STREQUAL "Fujitsu")))
    target_compile_features(kokkosprinter-tool PUBLIC cxx_std_14)
  endif()

  kokkos_add_test_executable(ProfilingAllCalls tools/TestAllCalls.cpp)

  kokkos_add_test_executable(ToolsInitialization UnitTestMain.cpp tools/TestToolsInitialization.cpp)

  set(ADDRESS_REGEX "0x[0-9a-f]*")
  set(MEMSPACE_REGEX "[HC][ou][sd][ta][a-zA-Z]*")
  set(SIZE_REGEX "[0-9]*")
  set(SKIP_SCRATCH_INITIALIZATION_REGEX ".*")

  # check help works via environment variable
  kokkos_add_test(
    SKIP_TRIBITS
    NAME
    ProfilingTestLibraryLoadHelp
    EXE
    ProfilingAllCalls
    TOOL
    kokkosprinter-tool
    ARGS
    --kokkos-tools-help
    PASS_REGULAR_EXPRESSION
    "kokkosp_init_library::kokkosp_print_help:Kokkos_ProfilingAllCalls::kokkosp_finalize_library::"
  )

  # check help works via direct library specification
  kokkos_add_test(
    SKIP_TRIBITS
    NAME
    ProfilingTestLibraryCmdLineHelp
    EXE
    ProfilingAllCalls
    ARGS
    --kokkos-tools-help
    --kokkos-tools-libs=$<TARGET_FILE:kokkosprinter-tool>
    PASS_REGULAR_EXPRESSION
    "kokkosp_init_library::kokkosp_print_help:Kokkos_ProfilingAllCalls::kokkosp_finalize_library::"
  )

  kokkos_add_test(
    SKIP_TRIBITS
    NAME
    ProfilingTestLibraryLoad
    EXE
    ProfilingAllCalls
    TOOL
    kokkosprinter-tool
    ARGS
    --kokkos-tools-args="-c test delimit"
    PASS_REGULAR_EXPRESSION
    "kokkosp_init_library::kokkosp_parse_args:4:Kokkos_ProfilingAllCalls:-c:test:delimit::.*::kokkosp_allocate_data:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::kokkosp_begin_parallel_for:Kokkos::View::initialization [[]source] via memset:[0-9]+:0::kokkosp_end_parallel_for:0::kokkosp_allocate_data:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:40::kokkosp_begin_parallel_for:Kokkos::View::initialization [[]destination] via memset:[0-9]+:0::kokkosp_end_parallel_for:0::kokkosp_begin_deep_copy:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::.*kokkosp_end_deep_copy::kokkosp_begin_parallel_for:parallel_for:${SIZE_REGEX}:0::kokkosp_end_parallel_for:0::kokkosp_begin_parallel_reduce:parallel_reduce:${SIZE_REGEX}:1${SKIP_SCRATCH_INITIALIZATION_REGEX}::kokkosp_end_parallel_reduce:1::kokkosp_begin_parallel_scan:parallel_scan:${SIZE_REGEX}:2::kokkosp_end_parallel_scan:2::kokkosp_push_profile_region:push_region::kokkosp_pop_profile_region::kokkosp_create_profile_section:created_section:3::kokkosp_start_profile_section:3::kokkosp_stop_profile_section:3::kokkosp_destroy_profile_section:3::kokkosp_profile_event:profiling_event::kokkosp_declare_metadata:dogs:good::kokkosp_deallocate_data:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:40::kokkosp_deallocate_data:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::kokkosp_finalize_library::"
  )

  # Above will test that leading/trailing quotes are stripped bc ctest cmd args is:
  #       "--kokkos-tools-args="-c test delimit""
  # The bracket argument syntax: [=[ and ]=] used below ensures it is treated as
  # a single argument:
  #       "--kokkos-tools-args=-c test delimit"
  #
  # https://cmake.org/cmake/help/latest/manual/cmake-language.7.html#bracket-argument
  #
  kokkos_add_test(
    SKIP_TRIBITS
    NAME
    ProfilingTestLibraryCmdLine
    EXE
    ProfilingAllCalls
    ARGS
    [=[--kokkos-tools-args=-c test delimit]=]
    --kokkos-tools-libs=$<TARGET_FILE:kokkosprinter-tool>
    PASS_REGULAR_EXPRESSION
    "kokkosp_init_library::kokkosp_parse_args:4:Kokkos_ProfilingAllCalls:-c:test:delimit::.*::kokkosp_allocate_data:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::kokkosp_begin_parallel_for:Kokkos::View::initialization [[]source] via memset:[0-9]+:0::kokkosp_end_parallel_for:0::kokkosp_allocate_data:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:40::kokkosp_begin_parallel_for:Kokkos::View::initialization [[]destination] via memset:[0-9]+:0::kokkosp_end_parallel_for:0::kokkosp_begin_deep_copy:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::.*kokkosp_end_deep_copy::kokkosp_begin_parallel_for:parallel_for:${SIZE_REGEX}:0::kokkosp_end_parallel_for:0::kokkosp_begin_parallel_reduce:parallel_reduce:${SIZE_REGEX}:1${SKIP_SCRATCH_INITIALIZATION_REGEX}::kokkosp_end_parallel_reduce:1::kokkosp_begin_parallel_scan:parallel_scan:${SIZE_REGEX}:2::kokkosp_end_parallel_scan:2::kokkosp_push_profile_region:push_region::kokkosp_pop_profile_region::kokkosp_create_profile_section:created_section:3::kokkosp_start_profile_section:3::kokkosp_stop_profile_section:3::kokkosp_destroy_profile_section:3::kokkosp_profile_event:profiling_event::kokkosp_declare_metadata:dogs:good::kokkosp_deallocate_data:${MEMSPACE_REGEX}:destination:${ADDRESS_REGEX}:40::kokkosp_deallocate_data:${MEMSPACE_REGEX}:source:${ADDRESS_REGEX}:40::kokkosp_finalize_library::"
  )
endif() #KOKKOS_ENABLE_LIBDL
kokkos_add_test_executable(
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
kokkos_set_exe_property(StackTraceTestExec ENABLE_EXPORTS ON)

kokkos_add_test(NAME CoreUnitTest_StackTraceTest EXE StackTraceTestExec FAIL_REGULAR_EXPRESSION "FAILED")

if(KOKKOS_ENABLE_HWLOC)
  kokkos_add_executable_and_test(CoreUnitTest_HWLOC SOURCES UnitTestMain.cpp TestHWLOC.cpp)
endif()

function(KOKKOS_ADD_INCREMENTAL_TEST DEVICE)
  kokkos_option(${DEVICE}_EXCLUDE_TESTS "" STRING "Incremental test exclude list")
  # Add unit test main
  set(${DEVICE}_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMainInit.cpp)

  # Iterate over incremental tests in directory

  append_glob(INCREMENTAL_FILE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/incremental/*.hpp)

  set(DEVICE_NAME ${KOKKOS_${DEVICE}_NAME})
  foreach(CURRENT_FILE_PATH ${INCREMENTAL_FILE_LIST})
    get_filename_component(CURRENT_FILE_NAME ${CURRENT_FILE_PATH} NAME)
    string(REPLACE ".hpp" "" CURRENT_TEST_NAME ${CURRENT_FILE_NAME})
    if(NOT CURRENT_TEST_NAME IN_LIST Kokkos_${DEVICE}_EXCLUDE_TESTS)
      set(CURRENT_TEST_OUTPUT_FILENAME ${CURRENT_TEST_NAME}_${DEVICE})
      file(STRINGS ${CURRENT_FILE_PATH} CURRENT_REQUIRED_FEATURE_LINE REGEX "Kokkos_Feature_Level_Required")
      # From each test get level implementation required
      string(REGEX REPLACE ".*Kokkos_Feature_Level_Required:" "" CURRENT_REQUIRED_FEATURE_LEVEL
                           ${CURRENT_REQUIRED_FEATURE_LINE}
      )
      # Cross-reference list of dependencies with selected feature list > matching feature test files are added to test applications
      if(KOKKOS_${DEVICE}_FEATURE_LEVEL GREATER_EQUAL CURRENT_REQUIRED_FEATURE_LEVEL)
        configure_file(
          IncrementalTest.cpp.in ${CMAKE_BINARY_DIR}/core/unit_test/generated/${CURRENT_TEST_OUTPUT_FILENAME}.cpp
        )
        set(${DEVICE}_SOURCES ${${DEVICE}_SOURCES};
                              ${CMAKE_BINARY_DIR}/core/unit_test/generated/${CURRENT_TEST_OUTPUT_FILENAME}.cpp
        )
      endif()
    endif()
  endforeach()

  string(TOUPPER ${DEVICE} UC_DEVICE)

  kokkos_option(ENABLE_${UC_DEVICE} ON BOOL "ENABLE ${UC_DEVICE}")

  kokkos_add_executable_and_test(IncrementalTest_${DEVICE} SOURCES ${${DEVICE}_SOURCES})

  set(EXE_NAME ${PACKAGE_NAME}_IncrementalTest_${DEVICE})
  # Check that the target was actually created because in a TribITS build
  # where only tests marked as PERFORMANCE enabled it would not be.
  if(TARGET ${EXE_NAME})
    target_include_directories(${EXE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/incremental)
  endif()

endfunction()

foreach(DEVICE ${KOKKOS_ENABLED_DEVICES})
  kokkos_add_incremental_test(${DEVICE})
endforeach()

kokkos_add_executable_and_test(CoreUnitTest_CTestDevice SOURCES UnitTestMain.cpp TestCTestDevice.cpp)

kokkos_add_executable_and_test(
  CoreUnitTest_CMakePassCmdLineArgs SOURCES UnitTest_CMakePassCmdLineArgs.cpp ARGS "one 2 THREE"
)

kokkos_add_executable(CoreUnitTest_CMakeTriBITSCompatibility SOURCES UnitTest_CMakeTriBITSCompatibility.cpp)

kokkos_add_test(NAME CoreUnitTest_CMakeTriBITSCompatibilityWillFail EXE CoreUnitTest_CMakeTriBITSCompatibility)
set_property(TEST Kokkos_CoreUnitTest_CMakeTriBITSCompatibilityWillFail PROPERTY WILL_FAIL TRUE)

set(Kokkos_CoreUnitTest_CMakeTriBITSCompatibilityDisable_DISABLE ON)
kokkos_add_test(NAME CoreUnitTest_CMakeTriBITSCompatibilityDisable EXE CoreUnitTest_CMakeTriBITSCompatibility)

set(Kokkos_CoreUnitTest_CMakeTriBITSCompatibilityExtraArgs_EXTRA_ARGS "--kokkos-test-tribits-compatibility=1")
kokkos_add_test(NAME CoreUnitTest_CMakeTriBITSCompatibilityExtraArgs EXE CoreUnitTest_CMakeTriBITSCompatibility)

set(Kokkos_CoreUnitTest_CMakeTriBITSCompatibilityEnvironment_ENVIRONMENT "KOKKOS_TEST_TRIBITS_COMPATIBILITY=1")
kokkos_add_test(NAME CoreUnitTest_CMakeTriBITSCompatibilityEnvironment EXE CoreUnitTest_CMakeTriBITSCompatibility)

set_source_files_properties(UnitTest_DeviceAndThreads.cpp PROPERTIES LANGUAGE ${KOKKOS_COMPILE_LANGUAGE})
add_executable(Kokkos_CoreUnitTest_DeviceAndThreads UnitTest_DeviceAndThreads.cpp)
target_link_libraries(Kokkos_CoreUnitTest_DeviceAndThreads Kokkos::kokkoscore)
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND AND Python3_VERSION VERSION_GREATER_EQUAL 3.7)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.20)
    set(USE_SOURCE_PERMISSIONS_WHEN_SUPPORTED USE_SOURCE_PERMISSIONS)
  endif()
  file(
    GENERATE
    OUTPUT $<TARGET_FILE_DIR:Kokkos_CoreUnitTest_DeviceAndThreads>/TestDeviceAndThreads.py
    INPUT TestDeviceAndThreads.py
    ${USE_SOURCE_PERMISSIONS_WHEN_SUPPORTED}
  )
  add_test(NAME Kokkos_CoreUnitTest_DeviceAndThreads
           COMMAND ${Python3_EXECUTABLE}
                   $<TARGET_FILE_DIR:Kokkos_CoreUnitTest_DeviceAndThreads>/TestDeviceAndThreads.py -v
  )
endif()

if(KOKKOS_ENABLE_HEADER_SELF_CONTAINMENT_TESTS AND NOT WIN32)
  add_subdirectory(headers_self_contained)
endif()
