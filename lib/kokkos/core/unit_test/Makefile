KOKKOS_PATH = ../..

GTEST_PATH = ../../tpls/gtest

vpath %.cpp ${KOKKOS_PATH}/core/unit_test
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/default
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/serial
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/threads
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/openmp
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/openmptarget
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/hip
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/hpx
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/cuda


TEST_HEADERS = $(wildcard $(KOKKOS_PATH)/core/unit_test/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOS_PATH)/core/unit_test/*/*.hpp)

default: build_all
	echo "End Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
  CXX = $(KOKKOS_PATH)/bin/nvcc_wrapper
else
  CXX = g++
endif

CXXFLAGS = -O3
LINK ?= $(CXX)
LDFLAGS ?=
override LDFLAGS += -lpthread

include $(KOKKOS_PATH)/Makefile.kokkos

KOKKOS_CXXFLAGS += -I$(GTEST_PATH) -I${KOKKOS_PATH}/core/unit_test -I${KOKKOS_PATH}/core/unit_test/category_files

TEST_TARGETS =
TARGETS =

KOKKOS_INTERNAL_HAS_OPTIMIZATIONS := $(call kokkos_has_string,$(KOKKOS_CXXFLAGS),O3)
ifneq ($(KOKKOS_INTERNAL_HAS_OPTIMIZATIONS), 1)
   KOKKOS_INTERNAL_HAS_OPTIMIZATIONS := $(call kokkos_has_string,$(KOKKOS_CXXFLAGS),O2)
   ifneq ($(KOKKOS_INTERNAL_HAS_OPTIMIZATIONS), 1)
      KOKKOS_INTERNAL_HAS_OPTIMIZATIONS := $(call kokkos_has_string,$(CXXFLAGS),O3)
      ifneq ($(KOKKOS_INTERNAL_HAS_OPTIMIZATIONS), 1)
         KOKKOS_INTERNAL_HAS_OPTIMIZATIONS := $(call kokkos_has_string,$(CXXFLAGS),O2)
      endif
   endif
endif
KOKKOS_INTERNAL_USE_RDYNAMIC := $(call kokkos_has_string,$(KOKKOS_CXXFLAGS),rdynamic)
ifneq ($(KOKKOS_INTERNAL_USE_RDYNAMIC), 1)
   KOKKOS_INTERNAL_USE_RDYNAMIC := $(call kokkos_has_string,$(CXXFLAGS),rdynamic)
endif

ifeq ($(KOKKOS_INTERNAL_USE_RDYNAMIC),1) 
   ifneq ($(KOKKOS_INTERNAL_HAS_OPTIMIZATIONS),1) 
      STACK_TRACE_TERMINATE_FILTER :=_dynamic
   else
      STACK_TRACE_TERMINATE_FILTER :=
   endif
else
   STACK_TRACE_TERMINATE_FILTER :=
endif

TESTS = AtomicOperations_int AtomicOperations_unsignedint AtomicOperations_longint AtomicOperations_unsignedlongint AtomicOperations_longlongint AtomicOperations_double AtomicOperations_float AtomicOperations_complexdouble AtomicOperations_complexfloat AtomicViews Atomics BlockSizeDeduction Concepts Complex Crs DeepCopyAlignment FunctorAnalysis Init LocalDeepCopy MDRange_a MDRange_b MDRange_c MDRange_d MDRange_e MDRange_f Other RangePolicy RangePolicyRequire Reductions Reducers_a Reducers_b Reducers_c Reducers_d Reducers_e Reductions_DeviceView Scan SharedAlloc TeamBasic TeamReductionScan TeamScratch TeamTeamSize TeamVectorRange UniqueToken ViewAPI_a ViewAPI_b ViewAPI_c ViewAPI_d ViewAPI_e ViewCopy_a ViewCopy_b ViewLayoutStrideAssignment ViewMapping_a ViewMapping_b ViewMapping_subview ViewOfClass WorkGraph View_64bit ViewResize

tmp := $(foreach device, $(KOKKOS_DEVICELIST), \
  tmp2 := $(foreach test, $(TESTS), \
    $(if $(filter Test$(device)_$(test).cpp, $(shell ls Test$(device)_$(test).cpp 2>/dev/null)),,\
      $(shell echo "\#include <Test"$(device)"_Category.hpp>" > Test$(device)_$(test).cpp); \
      $(shell echo "\#include <Test"$(test)".hpp>" >> Test$(device)_$(test).cpp); \
     ) \
  ) \
)

GPU_SPACE_TESTS = SharedAlloc ViewAPI_a ViewAPI_b ViewAPI_c ViewAPI_d ViewAPI_e ViewCopy_a ViewCopy_b ViewMapping_a ViewMapping_b ViewMapping_subview

SUBVIEW_TESTS = SubView_a SubView_b SubView_c01 SubView_c02 SubView_c03 SubView_c04 SubView_c05 SubView_c06 SubView_c07 SubView_c08 SubView_c09 SubView_c10 SubView_c11 SubView_c12 SubView_c13

KOKKOS_SUBVIEW_DEVICELIST := $(filter-out Cuda, $(KOKKOS_DEVICELIST))

tmp := $(foreach device, $(KOKKOS_SUBVIEW_DEVICELIST), \
  tmp2 := $(foreach test, $(SUBVIEW_TESTS), \
    $(if $(filter Test$(device)_$(test).cpp, $(shell ls Test$(device)_$(test).cpp 2>/dev/null)),, \
      $(shell echo "\#include <Test"$(device)"_Category.hpp>" > Test$(device)_$(test).cpp); \
      $(shell echo "\#include <Test"$(test)".hpp>" >> Test$(device)_$(test).cpp); \
     ) \
  )\
)

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
    tmp2 := $(foreach test, $(SUBVIEW_TESTS), \
      $(if $(filter TestCuda_$(test).cpp, $(shell ls TestCuda_$(test).cpp 2>/dev/null)),,\
        $(shell echo "\#include <TestCudaUVM_Category.hpp>" > TestCuda_$(test).cpp); \
        $(shell echo "\#include <Test"$(test)".hpp>" >> TestCuda_$(test).cpp); \
      )\
    )

    GPU_SPACES = CudaHostPinned CudaUVM
    tmp := $(foreach space, $(GPU_SPACES), \
      tmp2 := $(foreach test, $(GPU_SPACE_TESTS), \
        $(if $(filter Test$(space)_$(test).cpp, $(shell ls Test$(space)_$(test).cpp 2>/dev/null)),,\
          $(shell echo "\#include <Test$(space)_Category.hpp>" > Test$(space)_$(test).cpp); \
          $(shell echo "\#include <Test"$(test)".hpp>" >> Test$(space)_$(test).cpp); \
        )\
      )\
    )

    OBJ_CUDA = UnitTestMainInit.o gtest-all.o
    OBJ_CUDA += TestCuda_Init.o
    OBJ_CUDA += TestCuda_SharedAlloc.o TestCudaUVM_SharedAlloc.o TestCudaHostPinned_SharedAlloc.o
    OBJ_CUDA += TestCuda_RangePolicy.o TestCuda_RangePolicyRequire.o
    OBJ_CUDA += TestCuda_ViewAPI_a.o TestCuda_ViewAPI_b.o TestCuda_ViewAPI_c.o TestCuda_ViewAPI_d.o TestCuda_ViewAPI_e.o TestCuda_ViewCopy_a.o TestCuda_ViewCopy_b.o
    OBJ_CUDA += TestCuda_DeepCopyAlignment.o
    OBJ_CUDA += TestCuda_ViewMapping_a.o TestCuda_ViewMapping_b.o TestCuda_ViewMapping_subview.o TestCuda_ViewResize.o TestCuda_ViewLayoutStrideAssignment.o
    OBJ_CUDA += TestCudaUVM_ViewAPI_a.o TestCudaUVM_ViewAPI_b.o TestCudaUVM_ViewAPI_c.o TestCudaUVM_ViewAPI_d.o TestCudaUVM_ViewAPI_e.o
    OBJ_CUDA += TestCudaUVM_ViewCopy_a.o TestCudaUVM_ViewCopy_b.o
    OBJ_CUDA += TestCudaUVM_ViewMapping_a.o TestCudaUVM_ViewMapping_b.o TestCudaUVM_ViewMapping_subview.o
    OBJ_CUDA += TestCudaHostPinned_ViewAPI_a.o TestCudaHostPinned_ViewAPI_b.o TestCudaHostPinned_ViewAPI_c.o TestCudaHostPinned_ViewAPI_d.o TestCudaHostPinned_ViewAPI_e.o
    OBJ_CUDA += TestCudaHostPinned_ViewCopy_a.o TestCudaHostPinned_ViewCopy_b.o
    OBJ_CUDA += TestCudaHostPinned_ViewMapping_a.o TestCudaHostPinned_ViewMapping_b.o TestCudaHostPinned_ViewMapping_subview.o
    OBJ_CUDA += TestCuda_View_64bit.o
    OBJ_CUDA += TestCuda_ViewOfClass.o
    OBJ_CUDA += TestCuda_SubView_a.o TestCuda_SubView_b.o
    OBJ_CUDA += TestCuda_SubView_c01.o TestCuda_SubView_c02.o TestCuda_SubView_c03.o
    OBJ_CUDA += TestCuda_SubView_c04.o TestCuda_SubView_c05.o TestCuda_SubView_c06.o
    OBJ_CUDA += TestCuda_SubView_c07.o TestCuda_SubView_c08.o TestCuda_SubView_c09.o
    OBJ_CUDA += TestCuda_SubView_c10.o TestCuda_SubView_c11.o TestCuda_SubView_c12.o
    OBJ_CUDA += TestCuda_SubView_c13.o
    OBJ_CUDA += TestCuda_Reductions.o TestCuda_Scan.o
    OBJ_CUDA += TestCuda_Reductions_DeviceView.o
    OBJ_CUDA += TestCuda_Reducers_a.o TestCuda_Reducers_b.o TestCuda_Reducers_c.o TestCuda_Reducers_d.o TestCuda_Reducers_e.o
    OBJ_CUDA += TestCuda_Complex.o
    OBJ_CUDA += TestCuda_AtomicOperations_int.o TestCuda_AtomicOperations_unsignedint.o TestCuda_AtomicOperations_longint.o
    OBJ_CUDA += TestCuda_AtomicOperations_unsignedlongint.o TestCuda_AtomicOperations_longlongint.o TestCuda_AtomicOperations_double.o TestCuda_AtomicOperations_float.o
    OBJ_CUDA += TestCuda_AtomicOperations_complexfloat.o TestCuda_AtomicOperations_complexdouble.o
    OBJ_CUDA += TestCuda_AtomicViews.o TestCuda_Atomics.o
    OBJ_CUDA += TestCuda_TeamBasic.o TestCuda_TeamScratch.o
    OBJ_CUDA += TestCuda_TeamReductionScan.o TestCuda_TeamTeamSize.o
    OBJ_CUDA += TestCuda_TeamVectorRange.o
    OBJ_CUDA += TestCuda_Other.o
    OBJ_CUDA += TestCuda_MDRange_a.o TestCuda_MDRange_b.o TestCuda_MDRange_c.o TestCuda_MDRange_d.o TestCuda_MDRange_e.o
    OBJ_CUDA += TestCuda_Crs.o
    OBJ_CUDA += TestCuda_Task.o TestCuda_WorkGraph.o
    OBJ_CUDA += TestCuda_Spaces.o
    OBJ_CUDA += TestCuda_UniqueToken.o
    OBJ_CUDA += TestCuda_LocalDeepCopy.o
    OBJ_CUDA += TestCuda_DebugSerialExecution.o
    OBJ_CUDA += TestCuda_DebugPinUVMSpace.o
    OBJ_CUDA += TestCuda_TeamScratchStreams.o

    TARGETS += KokkosCore_UnitTest_Cuda
    TARGETS += KokkosCore_UnitTest_CudaInterOpInit
    TARGETS += KokkosCore_UnitTest_CudaInterOpStreams
    TEST_TARGETS += test-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_THREADS), 1)
    OBJ_THREADS = UnitTestMainInit.o gtest-all.o
    OBJ_THREADS += TestThreads_Init.o
    OBJ_THREADS += TestThreads_SharedAlloc.o
    OBJ_THREADS += TestThreads_RangePolicy.o TestThreads_RangePolicyRequire.o
    OBJ_THREADS += TestThreads_View_64bit.o
    OBJ_THREADS += TestThreads_ViewAPI_a.o TestThreads_ViewAPI_b.o TestThreads_ViewAPI_c.o TestThreads_ViewAPI_d.o TestThreads_ViewAPI_e.o
    OBJ_THREADS += TestThreads_ViewCopy_a.o TestThreads_ViewCopy_b.o
    OBJ_THREADS += TestThreads_DeepCopyAlignment.o
    OBJ_THREADS += TestThreads_ViewMapping_a.o TestThreads_ViewMapping_b.o TestThreads_ViewMapping_subview.o TestThreads_ViewResize.o TestThreads_ViewLayoutStrideAssignment.o
    OBJ_THREADS += TestThreads_ViewOfClass.o
    OBJ_THREADS += TestThreads_SubView_a.o TestThreads_SubView_b.o
    OBJ_THREADS += TestThreads_SubView_c01.o TestThreads_SubView_c02.o TestThreads_SubView_c03.o
    OBJ_THREADS += TestThreads_SubView_c04.o TestThreads_SubView_c05.o TestThreads_SubView_c06.o
    OBJ_THREADS += TestThreads_SubView_c07.o TestThreads_SubView_c08.o TestThreads_SubView_c09.o
    OBJ_THREADS += TestThreads_SubView_c10.o TestThreads_SubView_c11.o TestThreads_SubView_c12.o
    OBJ_THREADS += TestThreads_Reductions.o TestThreads_Scan.o
    OBJ_THREADS += TestThreads_Reductions_DeviceView.o
    OBJ_THREADS += TestThreads_Reducers_a.o TestThreads_Reducers_b.o TestThreads_Reducers_c.o TestThreads_Reducers_d.o TestThreads_Reducers_e.o
    OBJ_THREADS += TestThreads_Complex.o
    OBJ_THREADS += TestThreads_AtomicOperations_int.o TestThreads_AtomicOperations_unsignedint.o TestThreads_AtomicOperations_longint.o
    OBJ_THREADS += TestThreads_AtomicOperations_unsignedlongint.o TestThreads_AtomicOperations_longlongint.o TestThreads_AtomicOperations_double.o TestThreads_AtomicOperations_float.o
    OBJ_THREADS += TestThreads_AtomicOperations_complexfloat.o TestThreads_AtomicOperations_complexdouble.o
    OBJ_THREADS += TestThreads_AtomicViews.o TestThreads_Atomics.o
    OBJ_THREADS += TestThreads_TeamBasic.o TestThreads_TeamScratch.o TestThreads_TeamTeamSize.o
    OBJ_THREADS += TestThreads_TeamReductionScan.o
    OBJ_THREADS += TestThreads_TeamVectorRange.o
    OBJ_THREADS += TestThreads_Other.o
    OBJ_THREADS += TestThreads_MDRange_a.o TestThreads_MDRange_b.o TestThreads_MDRange_c.o TestThreads_MDRange_d.o TestThreads_MDRange_e.o
    OBJ_THREADS += TestThreads_LocalDeepCopy.o

    TARGETS += KokkosCore_UnitTest_Threads

    TEST_TARGETS += test-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
    OBJ_OPENMP = UnitTestMainInit.o gtest-all.o
    OBJ_OPENMP += TestOpenMP_Init.o
    OBJ_OPENMP += TestOpenMP_SharedAlloc.o
    OBJ_OPENMP += TestOpenMP_RangePolicy.o TestOpenMP_RangePolicyRequire.o
    OBJ_OPENMP += TestOpenMP_View_64bit.o
    OBJ_OPENMP += TestOpenMP_ViewAPI_a.o TestOpenMP_ViewAPI_b.o TestOpenMP_ViewAPI_c.o TestOpenMP_ViewAPI_d.o TestOpenMP_ViewAPI_e.o
    OBJ_OPENMP += TestOpenMP_DeepCopyAlignment.o TestOpenMP_ViewCopy_a.o TestOpenMP_ViewCopy_b.o
    OBJ_OPENMP += TestOpenMP_ViewMapping_a.o TestOpenMP_ViewMapping_b.o TestOpenMP_ViewMapping_subview.o TestOpenMP_ViewResize.o TestOpenMP_ViewLayoutStrideAssignment.o
    OBJ_OPENMP += TestOpenMP_ViewOfClass.o
    OBJ_OPENMP += TestOpenMP_SubView_a.o TestOpenMP_SubView_b.o
    OBJ_OPENMP += TestOpenMP_SubView_c01.o TestOpenMP_SubView_c02.o TestOpenMP_SubView_c03.o
    OBJ_OPENMP += TestOpenMP_SubView_c04.o TestOpenMP_SubView_c05.o TestOpenMP_SubView_c06.o
    OBJ_OPENMP += TestOpenMP_SubView_c07.o TestOpenMP_SubView_c08.o TestOpenMP_SubView_c09.o
    OBJ_OPENMP += TestOpenMP_SubView_c10.o TestOpenMP_SubView_c11.o TestOpenMP_SubView_c12.o
    OBJ_OPENMP += TestOpenMP_SubView_c13.o
    OBJ_OPENMP += TestOpenMP_Reductions.o TestOpenMP_Scan.o
    OBJ_OPENMP += TestOpenMP_Reductions_DeviceView.o
    OBJ_OPENMP += TestOpenMP_Reducers_a.o TestOpenMP_Reducers_b.o TestOpenMP_Reducers_c.o TestOpenMP_Reducers_d.o TestOpenMP_Reducers_e.o
    OBJ_OPENMP += TestOpenMP_Complex.o
    OBJ_OPENMP += TestOpenMP_AtomicOperations_int.o TestOpenMP_AtomicOperations_unsignedint.o TestOpenMP_AtomicOperations_longint.o
    OBJ_OPENMP += TestOpenMP_AtomicOperations_unsignedlongint.o TestOpenMP_AtomicOperations_longlongint.o TestOpenMP_AtomicOperations_double.o TestOpenMP_AtomicOperations_float.o
    OBJ_OPENMP += TestOpenMP_AtomicOperations_complexfloat.o TestOpenMP_AtomicOperations_complexdouble.o
    OBJ_OPENMP += TestOpenMP_AtomicViews.o TestOpenMP_Atomics.o
    OBJ_OPENMP += TestOpenMP_TeamBasic.o TestOpenMP_TeamScratch.o
    OBJ_OPENMP += TestOpenMP_TeamReductionScan.o TestOpenMP_TeamTeamSize.o
    OBJ_OPENMP += TestOpenMP_TeamVectorRange.o
    OBJ_OPENMP += TestOpenMP_Other.o
    OBJ_OPENMP += TestOpenMP_MDRange_a.o TestOpenMP_MDRange_b.o TestOpenMP_MDRange_c.o TestOpenMP_MDRange_d.o TestOpenMP_MDRange_e.o
    OBJ_OPENMP += TestOpenMP_Crs.o
    OBJ_OPENMP += TestOpenMP_Task.o TestOpenMP_WorkGraph.o
    OBJ_OPENMP += TestOpenMP_UniqueToken.o
    OBJ_OPENMP += TestOpenMP_LocalDeepCopy.o

    TARGETS += KokkosCore_UnitTest_OpenMP
    TARGETS += KokkosCore_UnitTest_OpenMPInterOp

    TEST_TARGETS += test-openmp
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMPTARGET), 1)
    OBJ_OPENMPTARGET = UnitTestMainInit.o gtest-all.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_Init.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_SharedAlloc.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_RangePolicy.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_ViewAPI_a.o TestOpenMPTarget_ViewAPI_b.o TestOpenMPTarget_ViewAPI_c.o TestOpenMPTarget_ViewAPI_d.o  #Some commented out code
    #OBJ_OPENMPTARGET += TestOpenMPTarget_ViewAPI_e.o TestOpenMPTarget_ViewCopy_a.o TestOpenMPTarget_ViewCopy_b.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_DeepCopyAlignment.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_ViewMapping_a.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_ViewMapping_b.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_ViewMapping_subview.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_ViewOfClass.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_a.o TestOpenMPTarget_SubView_b.o
    #The following subview tests need something like UVM:
    #OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_c01.o TestOpenMPTarget_SubView_c02.o TestOpenMPTarget_SubView_c03.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_c04.o TestOpenMPTarget_SubView_c05.o TestOpenMPTarget_SubView_c06.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_c07.o TestOpenMPTarget_SubView_c08.o TestOpenMPTarget_SubView_c09.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_c10.o TestOpenMPTarget_SubView_c11.o TestOpenMPTarget_SubView_c12.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_Reductions.o # Need custom reductions
    OBJ_OPENMPTARGET += TestOpenMPTarget_Reducers_a.o TestOpenMPTarget_Reducers_b.o TestOpenMPTarget_Reducers_c.o TestOpenMPTarget_Reducers_d.o TestOpenMPTarget_Reducers_e.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_Scan.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_Complex.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicOperations_int.o TestOpenMPTarget_AtomicOperations_unsignedint.o TestOpenMPTarget_AtomicOperations_longint.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicOperations_unsignedlongint.o TestOpenMPTarget_AtomicOperations_longlongint.o TestOpenMPTarget_AtomicOperations_double.o TestOpenMPTarget_AtomicOperations_float.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicOperations_complexfloat.o 
    #OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicOperations_complexdouble.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicViews.o
    OBJ_OPENMPTARGET += TestOpenMPTarget_Atomics.o # Commented Out Arbitrary Type Atomics
    #OBJ_OPENMPTARGET += TestOpenMPTarget_TeamBasic.o # There is still a static function in this
    #OBJ_OPENMPTARGET += TestOpenMPTarget_TeamScratch.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_TeamScan.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_TeamReductionScan.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_Other.o
    #OBJ_OPENMPTARGET += TestOpenMPTarget_MDRange_a.o TestOpenMPTarget_MDRange_b.o TestOpenMPTarget_MDRange_c.o TestOpenMPTarget_MDRange_d.o TestOpenMPTarget_MDRange_d.e
    #OBJ_OPENMPTARGET += TestOpenMPTarget_Task.o

    TARGETS += KokkosCore_UnitTest_OpenMPTarget

    TEST_TARGETS += test-openmptarget
endif

ifeq ($(KOKKOS_INTERNAL_USE_HIP), 1)
	GPU_SPACES = HIPHostPinned
	tmp := $(foreach space, $(GPU_SPACES), \
	  tmp2 := $(foreach test, $(GPU_SPACE_TESTS), \
	    $(if $(filter Test$(space)_$(test).cpp, $(shell ls Test$(space)_$(test).cpp 2>/dev/null)),,\
	      $(shell echo "\#include <Test$(space)_Category.hpp>" > Test$(space)_$(test).cpp); \
	      $(shell echo "\#include <Test"$(test)".hpp>" >> Test$(space)_$(test).cpp); \
	    )\
	  )\
	)

	OBJ_HIP = UnitTestMainInit.o gtest-all.o
	OBJ_HIP += TestHIP_Init.o
	OBJ_HIP += TestHIP_Reducers_a.o TestHIP_Reducers_b.o TestHIP_Reducers_c.o TestHIP_Reducers_d.o TestHIP_Reducers_e.o
	OBJ_HIP += TestHIP_Reductions.o
	OBJ_HIP += TestHIP_MDRange_a.o TestHIP_MDRange_b.o TestHIP_MDRange_c.o TestHIP_MDRange_d.o TestHIP_MDRange_e.o
	OBJ_HIP += TestHIP_Spaces.o
	OBJ_HIP += TestHIPHostPinned_ViewAPI_a.o TestHIPHostPinned_ViewAPI_b.o TestHIPHostPinned_ViewAPI_c.o TestHIPHostPinned_ViewAPI_d.o TestHIPHostPinned_ViewAPI_e.o
	OBJ_HIP += TestHIPHostPinned_ViewCopy_a.o TestHIPHostPinned_ViewCopy_b.o
	OBJ_HIP += TestHIPHostPinned_ViewMapping_a.o TestHIPHostPinned_ViewMapping_b.o TestHIPHostPinned_ViewMapping_subview.o

	TARGETS += KokkosCore_UnitTest_HIP

	TEST_TARGETS += test-hip
endif

ifeq ($(KOKKOS_INTERNAL_USE_HPX), 1)
	OBJ_HPX = UnitTestMainInit.o gtest-all.o
	OBJ_HPX += TestHPX_Init.o
	OBJ_HPX += TestHPX_SharedAlloc.o
	OBJ_HPX += TestHPX_RangePolicy.o TestHPX_RangePolicyRequire.o
	OBJ_HPX += TestHPX_View_64bit.o
	OBJ_HPX += TestHPX_ViewAPI_a.o TestHPX_ViewAPI_b.o TestHPX_ViewAPI_c.o TestHPX_ViewAPI_d.o TestHPX_ViewAPI_e.o
	OBJ_HPX += TestHPX_ViewCopy_a.o TestHPX_ViewCopy_b.o
	OBJ_HPX += TestHPX_ViewMapping_a.o TestHPX_ViewMapping_b.o TestHPX_ViewMapping_subview.o TestHPX_ViewResize.o
	OBJ_HPX += TestHPX_ViewOfClass.o
	OBJ_HPX += TestHPX_SubView_a.o TestHPX_SubView_b.o
	OBJ_HPX += TestHPX_SubView_c01.o TestHPX_SubView_c02.o TestHPX_SubView_c03.o
	OBJ_HPX += TestHPX_SubView_c04.o TestHPX_SubView_c05.o TestHPX_SubView_c06.o
	OBJ_HPX += TestHPX_SubView_c07.o TestHPX_SubView_c08.o TestHPX_SubView_c09.o
	OBJ_HPX += TestHPX_SubView_c10.o TestHPX_SubView_c11.o TestHPX_SubView_c12.o
	OBJ_HPX += TestHPX_SubView_c13.o
	OBJ_HPX += TestHPX_Reductions.o
	OBJ_HPX += TestHPX_Scan.o
	OBJ_HPX += TestHPX_Reducers_a.o TestHPX_Reducers_b.o TestHPX_Reducers_c.o TestHPX_Reducers_d.o TestHPX_Reducers_e.o
	OBJ_HPX += TestHPX_Complex.o
	OBJ_HPX += TestHPX_AtomicOperations_int.o TestHPX_AtomicOperations_unsignedint.o TestHPX_AtomicOperations_longint.o
	OBJ_HPX += TestHPX_AtomicOperations_unsignedlongint.o TestHPX_AtomicOperations_longlongint.o TestHPX_AtomicOperations_double.o TestHPX_AtomicOperations_float.o
	OBJ_HPX += TestHPX_AtomicViews.o TestHPX_Atomics.o
	OBJ_HPX += TestHPX_TeamBasic.o
	OBJ_HPX += TestHPX_TeamVectorRange.o
	OBJ_HPX += TestHPX_TeamScratch.o
	OBJ_HPX += TestHPX_TeamReductionScan.o
	OBJ_HPX += TestHPX_Other.o
	OBJ_HPX += TestHPX_MDRange_a.o TestHPX_MDRange_b.o TestHPX_MDRange_c.o TestHPX_MDRange_d.o TestHPX_MDRange_e.o
	OBJ_HPX += TestHPX_Crs.o
	OBJ_HPX += TestHPX_Task.o
	OBJ_HPX += TestHPX_WorkGraph.o
	OBJ_HPX += TestHPX_UniqueToken.o

	TARGETS += KokkosCore_UnitTest_HPX
  TARGETS += KokkosCore_UnitTest_HPXInterOp

	TEST_TARGETS += test-hpx
endif

ifeq ($(KOKKOS_INTERNAL_USE_SERIAL), 1)
    OBJ_SERIAL = UnitTestMainInit.o gtest-all.o
    OBJ_SERIAL += TestSerial_Init.o
    OBJ_SERIAL += TestSerial_SharedAlloc.o
    OBJ_SERIAL += TestSerial_RangePolicy.o TestSerial_RangePolicyRequire.o
    OBJ_SERIAL += TestSerial_View_64bit.o
    OBJ_SERIAL += TestSerial_ViewAPI_a.o TestSerial_ViewAPI_b.o TestSerial_ViewAPI_c.o TestSerial_ViewAPI_d.o TestSerial_ViewAPI_e.o
    OBJ_SERIAL += TestSerial_DeepCopyAlignment.o TestSerial_ViewCopy_a.o TestSerial_ViewCopy_b.o
    OBJ_SERIAL += TestSerial_ViewMapping_a.o TestSerial_ViewMapping_b.o TestSerial_ViewMapping_subview.o TestSerial_ViewResize.o TestSerial_ViewLayoutStrideAssignment.o
    OBJ_SERIAL += TestSerial_ViewOfClass.o
    OBJ_SERIAL += TestSerial_SubView_a.o TestSerial_SubView_b.o
    OBJ_SERIAL += TestSerial_SubView_c01.o TestSerial_SubView_c02.o TestSerial_SubView_c03.o
    OBJ_SERIAL += TestSerial_SubView_c04.o TestSerial_SubView_c05.o TestSerial_SubView_c06.o
    OBJ_SERIAL += TestSerial_SubView_c07.o TestSerial_SubView_c08.o TestSerial_SubView_c09.o
    OBJ_SERIAL += TestSerial_SubView_c10.o TestSerial_SubView_c11.o TestSerial_SubView_c12.o
    OBJ_SERIAL += TestSerial_SubView_c13.o
    OBJ_SERIAL += TestSerial_Reductions.o TestSerial_Scan.o
    OBJ_SERIAL += TestSerial_Reductions_DeviceView.o
    OBJ_SERIAL += TestSerial_Reducers_a.o TestSerial_Reducers_b.o TestSerial_Reducers_c.o TestSerial_Reducers_d.o TestSerial_Reducers_e.o
    OBJ_SERIAL += TestSerial_Complex.o
    OBJ_SERIAL += TestSerial_AtomicOperations_int.o TestSerial_AtomicOperations_unsignedint.o TestSerial_AtomicOperations_longint.o
    OBJ_SERIAL += TestSerial_AtomicOperations_unsignedlongint.o TestSerial_AtomicOperations_longlongint.o TestSerial_AtomicOperations_double.o TestSerial_AtomicOperations_float.o
    OBJ_SERIAL += TestSerial_AtomicOperations_complexfloat.o TestSerial_AtomicOperations_complexdouble.o
    OBJ_SERIAL += TestSerial_AtomicViews.o TestSerial_Atomics.o
    OBJ_SERIAL += TestSerial_TeamBasic.o TestSerial_TeamScratch.o
    OBJ_SERIAL += TestSerial_TeamVectorRange.o
    OBJ_SERIAL += TestSerial_TeamReductionScan.o TestSerial_TeamTeamSize.o
    OBJ_SERIAL += TestSerial_Other.o
    #HCC_WORKAROUND
    ifneq ($(KOKKOS_INTERNAL_COMPILER_HCC), 1)
        OBJ_SERIAL += TestSerial_MDRange_a.o TestSerial_MDRange_b.o TestSerial_MDRange_c.o TestSerial_MDRange_d.o TestSerial_MDRange_e.o
    endif
    OBJ_SERIAL += TestSerial_Crs.o
    OBJ_SERIAL += TestSerial_Task.o TestSerial_WorkGraph.o
    OBJ_SERIAL += TestSerial_LocalDeepCopy.o

    TARGETS += KokkosCore_UnitTest_Serial

    TEST_TARGETS += test-serial
endif

OBJ_HWLOC = TestHWLOC.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_HWLOC
TEST_TARGETS += test-hwloc

OBJ_DEFAULT = UnitTestMainInit.o gtest-all.o
ifneq ($(KOKKOS_INTERNAL_USE_OPENMPTARGET), 1)
ifneq ($(KOKKOS_INTERNAL_COMPILER_HCC), 1)
  OBJ_DEFAULT += TestDefaultDeviceType.o
  OBJ_DEFAULT += TestDefaultDeviceType_a1.o TestDefaultDeviceType_b1.o TestDefaultDeviceType_c1.o
  OBJ_DEFAULT += TestDefaultDeviceType_a2.o TestDefaultDeviceType_b2.o TestDefaultDeviceType_c2.o
  OBJ_DEFAULT += TestDefaultDeviceType_a3.o TestDefaultDeviceType_b3.o TestDefaultDeviceType_c3.o
  OBJ_DEFAULT += TestDefaultDeviceType_d.o
endif
endif

TARGETS += KokkosCore_UnitTest_Default
TEST_TARGETS += test-default

TARGETS += KokkosCore_UnitTest_PushFinalizeHook
TEST_TARGETS += test-push-finalize-hook

TARGETS += KokkosCore_UnitTest_PushFinalizeHook_terminate
TEST_TARGETS += test-push-finalize-hook-terminate

TARGETS += KokkosCore_UnitTest_StackTraceTestExec
TEST_TARGETS += test-stack-trace
TEST_TARGETS += test-stack-trace-terminate
TEST_TARGETS += test-stack-trace-generic-term

NUM_INITTESTS = 16
INITTESTS_NUMBERS := $(shell seq 1 ${NUM_INITTESTS})
INITTESTS_TARGETS := $(addprefix KokkosCore_UnitTest_DefaultDeviceTypeInit_,${INITTESTS_NUMBERS})
TARGETS += ${INITTESTS_TARGETS}
INITTESTS_TEST_TARGETS := $(addprefix test-default-init-,${INITTESTS_NUMBERS})
TEST_TARGETS += ${INITTESTS_TEST_TARGETS}

KokkosCore_UnitTest_Cuda: $(OBJ_CUDA) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_CUDA) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Cuda

KokkosCore_UnitTest_CudaInterOpInit: UnitTestMain.o gtest-all.o TestCuda_InterOp_Init.o $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) UnitTestMain.o gtest-all.o TestCuda_InterOp_Init.o $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_CudaInterOpInit
KokkosCore_UnitTest_CudaInterOpStreams: UnitTestMain.o gtest-all.o TestCuda_InterOp_Streams.o $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) UnitTestMain.o gtest-all.o TestCuda_InterOp_Streams.o $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_CudaInterOpStreams

KokkosCore_UnitTest_Threads: $(OBJ_THREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_THREADS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Threads

KokkosCore_UnitTest_OpenMP: $(OBJ_OPENMP) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_OPENMP) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_OpenMP

KokkosCore_UnitTest_OpenMPInterOp: UnitTestMain.o gtest-all.o TestOpenMP_InterOp.o $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) UnitTestMain.o gtest-all.o TestOpenMP_InterOp.o $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_OpenMPInterOp

KokkosCore_UnitTest_OpenMPTarget: $(OBJ_OPENMPTARGET) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_OPENMPTARGET) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_OpenMPTarget

KokkosCore_UnitTest_Serial: $(OBJ_SERIAL) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_SERIAL) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Serial

KokkosCore_UnitTest_HIP: $(OBJ_HIP) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_HIP) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_HIP

KokkosCore_UnitTest_HPX: $(OBJ_HPX) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_HPX) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_HPX

KokkosCore_UnitTest_HPXInterOp: UnitTestMain.o gtest-all.o TestHPX_InterOp.o $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) UnitTestMain.o gtest-all.o TestHPX_InterOp.o $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_HPXInterOp

KokkosCore_UnitTest_HWLOC: $(OBJ_HWLOC) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_HWLOC) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_HWLOC

KokkosCore_UnitTest_AllocationTracker: $(OBJ_ALLOCATIONTRACKER) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_ALLOCATIONTRACKER) $(KOKKOS_LIBS) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(LIB) -o KokkosCore_UnitTest_AllocationTracker

KokkosCore_UnitTest_Default: $(OBJ_DEFAULT) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_DEFAULT) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Default

KokkosCore_UnitTest_PushFinalizeHook: $(OBJ_DEFAULT) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_DEFAULT) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_PushFinalizeHook

KokkosCore_UnitTest_PushFinalizeHook_terminate: $(OBJ_DEFAULT) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_DEFAULT) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_PushFinalizeHook_terminate


${INITTESTS_TARGETS}: KokkosCore_UnitTest_DefaultDeviceTypeInit_%: TestDefaultDeviceTypeInit_%.o UnitTestMain.o gtest-all.o $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) TestDefaultDeviceTypeInit_$*.o UnitTestMain.o gtest-all.o $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_DefaultDeviceTypeInit_$*

KokkosCore_UnitTest_StackTraceTestExec: TestStackTrace.o  TestStackTrace_f0.o TestStackTrace_f1.o TestStackTrace_f2.o TestStackTrace_f3.o TestStackTrace_f4.o $(KOKKOS_LINK_DEPENDS) gtest-all.o
	$(LINK) $(EXTRA_PATH) TestStackTrace.o TestStackTrace_f0.o TestStackTrace_f1.o TestStackTrace_f2.o TestStackTrace_f3.o TestStackTrace_f4.o gtest-all.o $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_StackTraceTestExec

test-cuda: KokkosCore_UnitTest_Cuda
	./KokkosCore_UnitTest_Cuda
	./KokkosCore_UnitTest_CudaInterOpInit
	./KokkosCore_UnitTest_CudaInterOpStreams

test-threads: KokkosCore_UnitTest_Threads
	./KokkosCore_UnitTest_Threads

test-openmp: KokkosCore_UnitTest_OpenMP
	./KokkosCore_UnitTest_OpenMP
	./KokkosCore_UnitTest_OpenMPInterOp

test-openmptarget: KokkosCore_UnitTest_OpenMPTarget
	./KokkosCore_UnitTest_OpenMPTarget

test-serial: KokkosCore_UnitTest_Serial
	./KokkosCore_UnitTest_Serial

test-hip: KokkosCore_UnitTest_HIP
	./KokkosCore_UnitTest_HIP

test-hpx: KokkosCore_UnitTest_HPX
	./KokkosCore_UnitTest_HPX
	./KokkosCore_UnitTest_HPXInterOp

test-hwloc: KokkosCore_UnitTest_HWLOC
	./KokkosCore_UnitTest_HWLOC

test-allocationtracker: KokkosCore_UnitTest_AllocationTracker
	./KokkosCore_UnitTest_AllocationTracker

test-default: KokkosCore_UnitTest_Default
	./KokkosCore_UnitTest_Default

test-push-finalize-hook: KokkosCore_UnitTest_PushFinalizeHook
	./KokkosCore_UnitTest_PushFinalizeHook

test-push-finalize-hook-terminate: KokkosCore_UnitTest_PushFinalizeHook_terminate
	./KokkosCore_UnitTest_PushFinalizeHook_terminate
	
test-stack-trace: KokkosCore_UnitTest_StackTraceTestExec
	./KokkosCore_UnitTest_StackTraceTestExec --gtest_filter=*normal$(STACK_TRACE_TERMINATE_FILTER)

test-stack-trace-terminate: KokkosCore_UnitTest_StackTraceTestExec
	./KokkosCore_UnitTest_StackTraceTestExec --gtest_filter=*terminate$(STACK_TRACE_TERMINATE_FILTER)

test-stack-trace-generic-term: KokkosCore_UnitTest_StackTraceTestExec
	./KokkosCore_UnitTest_StackTraceTestExec --gtest_filter=*generic_term$(STACK_TRACE_TERMINATE_FILTER)


${INITTESTS_TEST_TARGETS}: test-default-init-%: KokkosCore_UnitTest_DefaultDeviceTypeInit_%
	./KokkosCore_UnitTest_DefaultDeviceTypeInit_$*

build_all: $(TARGETS)

test: $(TEST_TARGETS)

clean: kokkos-clean
	rm -f *.o $(TARGETS) TestCuda*.cpp TestThreads*.cpp TestOpenMP*.cpp TestSerial*.cpp TestHIP*.cpp \
        TestOpenMPTarget*.cpp TestHPX*.cpp

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(TEST_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

gtest-all.o:$(GTEST_PATH)/gtest/gtest-all.cc
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $(GTEST_PATH)/gtest/gtest-all.cc
