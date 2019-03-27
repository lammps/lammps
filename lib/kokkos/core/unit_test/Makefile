KOKKOS_PATH = ../..

GTEST_PATH = ../../tpls/gtest

vpath %.cpp ${KOKKOS_PATH}/core/unit_test
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/default
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/serial
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/threads
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/openmp
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/openmptarget
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/qthreads
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/cuda
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/rocm


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

KOKKOS_CXXFLAGS += -I$(GTEST_PATH) -I${KOKKOS_PATH}/core/unit_test

TEST_TARGETS =
TARGETS =

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
	OBJ_CUDA = UnitTestMainInit.o gtest-all.o
	OBJ_CUDA += TestCuda_Init.o
	OBJ_CUDA += TestCuda_SharedAlloc.o TestCudaUVM_SharedAlloc.o TestCudaHostPinned_SharedAlloc.o
	OBJ_CUDA += TestCuda_RangePolicy.o
	OBJ_CUDA += TestCuda_ViewAPI_a.o TestCuda_ViewAPI_b.o TestCuda_ViewAPI_c.o TestCuda_ViewAPI_d.o TestCuda_ViewAPI_e.o
	OBJ_CUDA += TestCuda_ViewMapping_a.o TestCuda_ViewMapping_b.o TestCuda_ViewMapping_subview.o TestCuda_ViewLayoutStrideAssignment.o
	OBJ_CUDA += TestCudaUVM_ViewCopy.o TestCudaUVM_ViewAPI_a.o TestCudaUVM_ViewAPI_b.o TestCudaUVM_ViewAPI_c.o TestCudaUVM_ViewAPI_d.o TestCudaUVM_ViewAPI_e.o
	OBJ_CUDA += TestCudaUVM_ViewMapping_a.o TestCudaUVM_ViewMapping_b.o TestCudaUVM_ViewMapping_subview.o
	OBJ_CUDA += TestCudaHostPinned_ViewCopy.o TestCudaHostPinned_ViewAPI_a.o TestCudaHostPinned_ViewAPI_b.o TestCudaHostPinned_ViewAPI_c.o TestCudaHostPinned_ViewAPI_d.o TestCudaHostPinned_ViewAPI_e.o
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
	OBJ_CUDA += TestCuda_Reducers_a.o TestCuda_Reducers_b.o TestCuda_Reducers_c.o TestCuda_Reducers_d.o
	OBJ_CUDA += TestCuda_Complex.o
	OBJ_CUDA += TestCuda_AtomicOperations_int.o TestCuda_AtomicOperations_unsignedint.o TestCuda_AtomicOperations_longint.o 
	OBJ_CUDA += TestCuda_AtomicOperations_unsignedlongint.o TestCuda_AtomicOperations_longlongint.o TestCuda_AtomicOperations_double.o TestCuda_AtomicOperations_float.o
	OBJ_CUDA += TestCuda_AtomicViews.o TestCuda_Atomics.o
	OBJ_CUDA += TestCuda_Team.o TestCuda_TeamScratch.o
	OBJ_CUDA += TestCuda_TeamReductionScan.o TestCuda_TeamTeamSize.o
	OBJ_CUDA += TestCuda_Other.o
	OBJ_CUDA += TestCuda_MDRange_a.o TestCuda_MDRange_b.o TestCuda_MDRange_c.o TestCuda_MDRange_d.o TestCuda_MDRange_e.o
	OBJ_CUDA += TestCuda_Crs.o
	OBJ_CUDA += TestCuda_Task.o TestCuda_WorkGraph.o
	OBJ_CUDA += TestCuda_Spaces.o
	OBJ_CUDA += TestCuda_UniqueToken.o
	
	TARGETS += KokkosCore_UnitTest_Cuda
    TARGETS += KokkosCore_UnitTest_CudaInterOp
	TEST_TARGETS += test-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_ROCM), 1)
        OBJ_ROCM = UnitTestMainInit.o gtest-all.o
        OBJ_ROCM += TestROCm_Init.o
        OBJ_ROCM += TestROCm_Complex.o
        OBJ_ROCM += TestROCm_RangePolicy.o
	OBJ_ROCM += TestROCm_AtomicOperations_int.o TestROCm_AtomicOperations_unsignedint.o TestROCm_AtomicOperations_longint.o 
	OBJ_ROCM += TestROCm_AtomicOperations_unsignedlongint.o TestROCm_AtomicOperations_longlongint.o TestROCm_AtomicOperations_double.o TestROCm_AtomicOperations_float.o
        OBJ_ROCM += TestROCm_Atomics.o
        OBJ_ROCM += TestROCm_AtomicViews.o
        OBJ_ROCM += TestROCm_Other.o
        OBJ_ROCM += TestROCm_MDRange_a.o TestROCm_MDRange_b.o TestROCm_MDRange_c.o TestROCm_MDRange_d.o TestROCm_MDRange_e.o
        OBJ_ROCM += TestROCm_MDRangeReduce_a.o TestROCm_MDRangeReduce_b.o TestROCm_MDRangeReduce_c.o TestROCm_MDRangeReduce_d.o TestROCm_MDRangeReduce_e.o
	OBJ_ROCM += TestROCm_Reductions.o
	OBJ_ROCM += TestROCm_Reducers_a.o TestROCm_Reducers_b.o TestROCm_Reducers_c.o TestROCm_Reducers_d.o
        OBJ_ROCM += TestROCm_Scan.o
        OBJ_ROCM += TestROCm_SharedAlloc.o
        OBJ_ROCM += TestROCm_SubView_a.o
        OBJ_ROCM += TestROCm_SubView_b.o
        OBJ_ROCM += TestROCm_SubView_c01.o
        OBJ_ROCM += TestROCm_SubView_c02.o
        OBJ_ROCM += TestROCm_SubView_c03.o
        OBJ_ROCM += TestROCm_SubView_c04.o
        OBJ_ROCM += TestROCm_SubView_c05.o
        OBJ_ROCM += TestROCm_SubView_c06.o
        OBJ_ROCM += TestROCm_SubView_c07.o
        OBJ_ROCM += TestROCm_SubView_c08.o
        OBJ_ROCM += TestROCm_SubView_c09.o
        OBJ_ROCM += TestROCm_SubView_c10.o
        OBJ_ROCM += TestROCm_SubView_c11.o
        OBJ_ROCM += TestROCm_SubView_c12.o
        OBJ_ROCM += TestROCm_SubView_c13.o
        OBJ_ROCM += TestROCm_Team.o
        OBJ_ROCM += TestROCm_TeamReductionScan.o
        OBJ_ROCM += TestROCm_TeamScratch.o TestROCm_TeamTeamSize.o
        OBJ_ROCM += TestROCm_ViewAPI_a.o TestROCm_ViewAPI_b.o TestROCm_ViewAPI_c.o TestROCm_ViewAPI_d.o TestROCm_ViewAPI_e.o
        OBJ_ROCM += TestROCm_ViewMapping_a.o
        OBJ_ROCM += TestROCm_ViewMapping_b.o
        OBJ_ROCM += TestROCm_ViewMapping_subview.o
        OBJ_ROCM += TestROCmHostPinned_ViewCopy.o TestROCmHostPinned_ViewAPI_a.o TestROCmHostPinned_ViewAPI_b.o TestROCmHostPinned_ViewAPI_c.o TestROCmHostPinned_ViewAPI_d.o TestROCmHostPinned_ViewAPI_e.o
        OBJ_ROCM += TestROCmHostPinned_View_64bit.o
        OBJ_ROCM += TestROCmHostPinned_ViewMapping_a.o 
        OBJ_ROCM += TestROCmHostPinned_ViewMapping_b.o 
        OBJ_ROCM += TestROCmHostPinned_ViewMapping_subview.o
        OBJ_ROCM += TestROCm_ViewOfClass.o
        OBJ_ROCM += TestROCm_Spaces.o
        OBJ_ROCM += TestROCm_Crs.o
     
        TARGETS += KokkosCore_UnitTest_ROCm
        TEST_TARGETS += test-rocm
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)

	OBJ_THREADS = UnitTestMainInit.o gtest-all.o
	OBJ_THREADS += TestThreads_Init.o
	OBJ_THREADS += TestThreads_SharedAlloc.o
	OBJ_THREADS += TestThreads_RangePolicy.o
    OBJ_THREADS += TestThreads_View_64bit.o
	OBJ_THREADS += TestThreads_ViewAPI_a.o TestThreads_ViewAPI_b.o TestThreads_ViewAPI_c.o TestThreads_ViewAPI_d.o TestThreads_ViewAPI_e.o
	OBJ_THREADS += TestThreads_ViewMapping_a.o TestThreads_ViewMapping_b.o TestThreads_ViewMapping_subview.o TestThreads_ViewLayoutStrideAssignment.o
	OBJ_THREADS += TestThreads_ViewOfClass.o
	OBJ_THREADS += TestThreads_SubView_a.o TestThreads_SubView_b.o
	OBJ_THREADS += TestThreads_SubView_c01.o TestThreads_SubView_c02.o TestThreads_SubView_c03.o
	OBJ_THREADS += TestThreads_SubView_c04.o TestThreads_SubView_c05.o TestThreads_SubView_c06.o
	OBJ_THREADS += TestThreads_SubView_c07.o TestThreads_SubView_c08.o TestThreads_SubView_c09.o
	OBJ_THREADS += TestThreads_SubView_c10.o TestThreads_SubView_c11.o TestThreads_SubView_c12.o
	OBJ_THREADS += TestThreads_Reductions.o TestThreads_Scan.o
	OBJ_THREADS += TestThreads_Reductions_DeviceView.o
	OBJ_THREADS += TestThreads_Reducers_a.o TestThreads_Reducers_b.o TestThreads_Reducers_c.o TestThreads_Reducers_d.o
	OBJ_THREADS += TestThreads_Complex.o
	OBJ_THREADS += TestThreads_AtomicOperations_int.o TestThreads_AtomicOperations_unsignedint.o TestThreads_AtomicOperations_longint.o 
	OBJ_THREADS += TestThreads_AtomicOperations_unsignedlongint.o TestThreads_AtomicOperations_longlongint.o TestThreads_AtomicOperations_double.o TestThreads_AtomicOperations_float.o
	OBJ_THREADS += TestThreads_AtomicViews.o TestThreads_Atomics.o
	OBJ_THREADS += TestThreads_Team.o TestThreads_TeamScratch.o TestThreads_TeamTeamSize.o
	OBJ_THREADS += TestThreads_TeamReductionScan.o
	OBJ_THREADS += TestThreads_Other.o
	OBJ_THREADS += TestThreads_MDRange_a.o TestThreads_MDRange_b.o TestThreads_MDRange_c.o TestThreads_MDRange_d.o TestThreads_MDRange_e.o

	TARGETS += KokkosCore_UnitTest_Threads

	TEST_TARGETS += test-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
	OBJ_OPENMP = UnitTestMainInit.o gtest-all.o
	OBJ_OPENMP += TestOpenMP_Init.o
	OBJ_OPENMP += TestOpenMP_SharedAlloc.o
	OBJ_OPENMP += TestOpenMP_RangePolicy.o
    OBJ_OPENMP += TestOpenMP_View_64bit.o
	OBJ_OPENMP += TestOpenMP_ViewAPI_a.o TestOpenMP_ViewAPI_b.o TestOpenMP_ViewAPI_c.o TestOpenMP_ViewAPI_d.o TestOpenMP_ViewAPI_e.o
	OBJ_OPENMP += TestOpenMP_ViewMapping_a.o TestOpenMP_ViewMapping_b.o TestOpenMP_ViewMapping_subview.o TestOpenMP_ViewLayoutStrideAssignment.o
	OBJ_OPENMP += TestOpenMP_ViewOfClass.o
	OBJ_OPENMP += TestOpenMP_SubView_a.o TestOpenMP_SubView_b.o
	OBJ_OPENMP += TestOpenMP_SubView_c01.o TestOpenMP_SubView_c02.o TestOpenMP_SubView_c03.o
	OBJ_OPENMP += TestOpenMP_SubView_c04.o TestOpenMP_SubView_c05.o TestOpenMP_SubView_c06.o
	OBJ_OPENMP += TestOpenMP_SubView_c07.o TestOpenMP_SubView_c08.o TestOpenMP_SubView_c09.o
	OBJ_OPENMP += TestOpenMP_SubView_c10.o TestOpenMP_SubView_c11.o TestOpenMP_SubView_c12.o
	OBJ_OPENMP += TestOpenMP_SubView_c13.o
	OBJ_OPENMP += TestOpenMP_Reductions.o TestOpenMP_Scan.o
	OBJ_OPENMP += TestOpenMP_Reductions_DeviceView.o
	OBJ_OPENMP += TestOpenMP_Reducers_a.o TestOpenMP_Reducers_b.o TestOpenMP_Reducers_c.o TestOpenMP_Reducers_d.o
	OBJ_OPENMP += TestOpenMP_Complex.o
	OBJ_OPENMP += TestOpenMP_AtomicOperations_int.o TestOpenMP_AtomicOperations_unsignedint.o TestOpenMP_AtomicOperations_longint.o 
	OBJ_OPENMP += TestOpenMP_AtomicOperations_unsignedlongint.o TestOpenMP_AtomicOperations_longlongint.o TestOpenMP_AtomicOperations_double.o TestOpenMP_AtomicOperations_float.o
	OBJ_OPENMP += TestOpenMP_AtomicViews.o TestOpenMP_Atomics.o
	OBJ_OPENMP += TestOpenMP_Team.o TestOpenMP_TeamScratch.o
	OBJ_OPENMP += TestOpenMP_TeamReductionScan.o TestOpenMP_TeamTeamSize.o
	OBJ_OPENMP += TestOpenMP_Other.o
	OBJ_OPENMP += TestOpenMP_MDRange_a.o TestOpenMP_MDRange_b.o TestOpenMP_MDRange_c.o TestOpenMP_MDRange_d.o TestOpenMP_MDRange_e.o
	OBJ_OPENMP += TestOpenMP_Crs.o
	OBJ_OPENMP += TestOpenMP_Task.o TestOpenMP_WorkGraph.o
	OBJ_OPENMP += TestOpenMP_UniqueToken.o
	
	TARGETS += KokkosCore_UnitTest_OpenMP
    TARGETS += KokkosCore_UnitTest_OpenMPInterOp

	TEST_TARGETS += test-openmp
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMPTARGET), 1)
	OBJ_OPENMPTARGET = UnitTestMainInit.o gtest-all.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_Init.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_SharedAlloc.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_RangePolicy.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_ViewAPI_a.o TestOpenMPTarget_ViewAPI_b.o TestOpenMPTarget_ViewAPI_c.o TestOpenMPTarget_ViewAPI_d.o TestOpenMPTarget_ViewAPI_e.o #Some commented out code
	OBJ_OPENMPTARGET += TestOpenMPTarget_ViewMapping_a.o 
    OBJ_OPENMPTARGET += TestOpenMPTarget_ViewMapping_b.o 
    OBJ_OPENMPTARGET += TestOpenMPTarget_ViewMapping_subview.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_ViewOfClass.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_a.o TestOpenMPTarget_SubView_b.o
	#The following subview tests need something like UVM:
        #OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_c01.o TestOpenMPTarget_SubView_c02.o TestOpenMPTarget_SubView_c03.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_c04.o TestOpenMPTarget_SubView_c05.o TestOpenMPTarget_SubView_c06.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_c07.o TestOpenMPTarget_SubView_c08.o TestOpenMPTarget_SubView_c09.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_SubView_c10.o TestOpenMPTarget_SubView_c11.o TestOpenMPTarget_SubView_c12.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_Reductions.o # Need custom reductions
	#OBJ_OPENMPTARGET += TestOpenMPTarget_Reducers_a.o TestOpenMPTarget_Reducers_b.o TestOpenMPTarget_Reducers_c.o TestOpenMPTarget_Reducers_d.o
        #OBJ_OPENMPTARGET += TestOpenMPTarget_Scan.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_Complex.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicOperations_int.o TestOpenMPTarget_AtomicOperations_unsignedint.o TestOpenMPTarget_AtomicOperations_longint.o 
	OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicOperations_unsignedlongint.o TestOpenMPTarget_AtomicOperations_longlongint.o TestOpenMPTarget_AtomicOperations_double.o TestOpenMPTarget_AtomicOperations_float.o
        OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicViews.o 
        OBJ_OPENMPTARGET += TestOpenMPTarget_Atomics.o # Commented Out Arbitrary Type Atomics
	#OBJ_OPENMPTARGET += TestOpenMPTarget_Team.o # There is still a static function in this
        #OBJ_OPENMPTARGET += TestOpenMPTarget_TeamScratch.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_TeamReductionScan.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_Other.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_MDRange_a.o TestOpenMPTarget_MDRange_b.o TestOpenMPTarget_MDRange_c.o TestOpenMPTarget_MDRange_d.o TestOpenMPTarget_MDRange_d.e
	#OBJ_OPENMPTARGET += TestOpenMPTarget_Task.o
	
	TARGETS += KokkosCore_UnitTest_OpenMPTarget

	TEST_TARGETS += test-openmptarget

endif

ifeq ($(KOKKOS_INTERNAL_USE_QTHREADS), 1)
	OBJ_QTHREADS = TestQthreads_Other.o TestQthreads_Reductions.o TestQthreads_Atomics.o TestQthreads_Team.o
	OBJ_QTHREADS += TestQthreads_SubView_a.o TestQthreads_SubView_b.o
	OBJ_QTHREADS += TestQthreads_SubView_c01.o TestQthreads_SubView_c02.o TestQthreads_SubView_c03.o
	OBJ_QTHREADS += TestQthreads_SubView_c04.o TestQthreads_SubView_c05.o TestQthreads_SubView_c06.o
	OBJ_QTHREADS += TestQthreads_SubView_c07.o TestQthreads_SubView_c08.o TestQthreads_SubView_c09.o
	OBJ_QTHREADS += TestQthreads_SubView_c10.o TestQthreads_SubView_c11.o TestQthreads_SubView_c12.o
	OBJ_QTHREADS += TestQthreads_ViewAPI_a.o TestQthreads_ViewAPI_b.o TestQthreads_ViewAPI_c.o TestQthreads_ViewAPI_d.o TestQthreads_ViewAPI_e.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_Qthreads

	OBJ_QTHREADS2 = UnitTestMainInit.o gtest-all.o
	OBJ_QTHREADS2 += TestQthreads_Complex.o
	TARGETS += KokkosCore_UnitTest_Qthreads2

	TEST_TARGETS += test-qthreads
endif

ifeq ($(KOKKOS_INTERNAL_USE_SERIAL), 1)
        OBJ_SERIAL = UnitTestMainInit.o gtest-all.o
        OBJ_SERIAL += TestSerial_Init.o
        OBJ_SERIAL += TestSerial_SharedAlloc.o
        OBJ_SERIAL += TestSerial_RangePolicy.o
        OBJ_SERIAL += TestSerial_View_64bit.o
        OBJ_SERIAL += TestSerial_ViewAPI_a.o TestSerial_ViewAPI_b.o TestSerial_ViewAPI_c.o TestSerial_ViewAPI_d.o TestSerial_ViewAPI_e.o
        OBJ_SERIAL += TestSerial_ViewMapping_a.o TestSerial_ViewMapping_b.o TestSerial_ViewMapping_subview.o TestSerial_ViewLayoutStrideAssignment.o
        OBJ_SERIAL += TestSerial_ViewOfClass.o
        OBJ_SERIAL += TestSerial_SubView_a.o TestSerial_SubView_b.o
        OBJ_SERIAL += TestSerial_SubView_c01.o TestSerial_SubView_c02.o TestSerial_SubView_c03.o
        OBJ_SERIAL += TestSerial_SubView_c04.o TestSerial_SubView_c05.o TestSerial_SubView_c06.o
        OBJ_SERIAL += TestSerial_SubView_c07.o TestSerial_SubView_c08.o TestSerial_SubView_c09.o
        OBJ_SERIAL += TestSerial_SubView_c10.o TestSerial_SubView_c11.o TestSerial_SubView_c12.o
        OBJ_SERIAL += TestSerial_SubView_c13.o
        OBJ_SERIAL += TestSerial_Reductions.o TestSerial_Scan.o
	OBJ_SERIAL += TestSerial_Reductions_DeviceView.o
    	OBJ_SERIAL += TestSerial_Reducers_a.o TestSerial_Reducers_b.o TestSerial_Reducers_c.o TestSerial_Reducers_d.o
        OBJ_SERIAL += TestSerial_Complex.o
	    OBJ_SERIAL += TestSerial_AtomicOperations_int.o TestSerial_AtomicOperations_unsignedint.o TestSerial_AtomicOperations_longint.o 
	    OBJ_SERIAL += TestSerial_AtomicOperations_unsignedlongint.o TestSerial_AtomicOperations_longlongint.o TestSerial_AtomicOperations_double.o TestSerial_AtomicOperations_float.o
	    OBJ_SERIAL += TestSerial_AtomicViews.o TestSerial_Atomics.o
        OBJ_SERIAL += TestSerial_Team.o TestSerial_TeamScratch.o
        OBJ_SERIAL += TestSerial_TeamReductionScan.o TestSerial_TeamTeamSize.o
        OBJ_SERIAL += TestSerial_Other.o
        #HCC_WORKAROUND
        ifneq ($(KOKKOS_INTERNAL_COMPILER_HCC), 1)
        OBJ_SERIAL += TestSerial_MDRange_a.o TestSerial_MDRange_b.o TestSerial_MDRange_c.o TestSerial_MDRange_d.o TestSerial_MDRange_e.o
        endif
        OBJ_SERIAL += TestSerial_Crs.o
        OBJ_SERIAL += TestSerial_Task.o TestSerial_WorkGraph.o
	
	TARGETS += KokkosCore_UnitTest_Serial

	TEST_TARGETS += test-serial
endif

OBJ_HWLOC = TestHWLOC.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_HWLOC
TEST_TARGETS += test-hwloc

OBJ_HOST_BARRIER = TestHostBarrier.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_HostBarrier
TEST_TARGETS += test-host-barrier

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

NUM_INITTESTS = 16
INITTESTS_NUMBERS := $(shell seq 1 ${NUM_INITTESTS})
INITTESTS_TARGETS := $(addprefix KokkosCore_UnitTest_DefaultDeviceTypeInit_,${INITTESTS_NUMBERS})
TARGETS += ${INITTESTS_TARGETS}
INITTESTS_TEST_TARGETS := $(addprefix test-default-init-,${INITTESTS_NUMBERS})
TEST_TARGETS += ${INITTESTS_TEST_TARGETS}

KokkosCore_UnitTest_Cuda: $(OBJ_CUDA) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_CUDA) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Cuda
	
KokkosCore_UnitTest_CudaInterOp: UnitTestMain.o gtest-all.o TestCuda_InterOp.o
	$(LINK) $(EXTRA_PATH) UnitTestMain.o gtest-all.o TestCuda_InterOp.o $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_CudaInterOp
	
KokkosCore_UnitTest_ROCm: $(OBJ_ROCM) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_ROCM) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_ROCm

KokkosCore_UnitTest_Threads: $(OBJ_THREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_THREADS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Threads

KokkosCore_UnitTest_OpenMP: $(OBJ_OPENMP) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_OPENMP) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_OpenMP

KokkosCore_UnitTest_OpenMPInterOp: UnitTestMain.o gtest-all.o TestOpenMP_InterOp.o
	$(LINK) $(EXTRA_PATH) UnitTestMain.o gtest-all.o TestOpenMP_InterOp.o $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_OpenMPInterOp

KokkosCore_UnitTest_OpenMPTarget: $(OBJ_OPENMPTARGET) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_OPENMPTARGET) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_OpenMPTarget

KokkosCore_UnitTest_Serial: $(OBJ_SERIAL) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_SERIAL) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Serial

KokkosCore_UnitTest_Qthreads: $(OBJ_QTHREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_QTHREADS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Qthreads

KokkosCore_UnitTest_Qthreads2: $(OBJ_QTHREADS2) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_QTHREADS2) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Qthreads2

KokkosCore_UnitTest_HWLOC: $(OBJ_HWLOC) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_HWLOC) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_HWLOC

KokkosCore_UnitTest_HostBarrier: $(OBJ_HOST_BARRIER) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_HOST_BARRIER) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_HostBarrier

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

test-cuda: KokkosCore_UnitTest_Cuda
	./KokkosCore_UnitTest_Cuda
	./KokkosCore_UnitTest_CudaInterOp

test-rocm: KokkosCore_UnitTest_ROCm
	./KokkosCore_UnitTest_ROCm

test-threads: KokkosCore_UnitTest_Threads
	./KokkosCore_UnitTest_Threads

test-openmp: KokkosCore_UnitTest_OpenMP
	./KokkosCore_UnitTest_OpenMP
	./KokkosCore_UnitTest_OpenMPInterOp

test-openmptarget: KokkosCore_UnitTest_OpenMPTarget
	./KokkosCore_UnitTest_OpenMPTarget

test-serial: KokkosCore_UnitTest_Serial
	./KokkosCore_UnitTest_Serial

test-qthreads: KokkosCore_UnitTest_Qthreads KokkosCore_UnitTest_Qthreads2
	./KokkosCore_UnitTest_Qthreads
	./KokkosCore_UnitTest_Qthreads2

test-hwloc: KokkosCore_UnitTest_HWLOC
	./KokkosCore_UnitTest_HWLOC

test-host-barrier: KokkosCore_UnitTest_HostBarrier
	./KokkosCore_UnitTest_HostBarrier

test-allocationtracker: KokkosCore_UnitTest_AllocationTracker
	./KokkosCore_UnitTest_AllocationTracker

test-default: KokkosCore_UnitTest_Default
	./KokkosCore_UnitTest_Default

test-push-finalize-hook: KokkosCore_UnitTest_PushFinalizeHook
	./KokkosCore_UnitTest_PushFinalizeHook

test-push-finalize-hook-terminate: KokkosCore_UnitTest_PushFinalizeHook_terminate
	./KokkosCore_UnitTest_PushFinalizeHook_terminate

${INITTESTS_TEST_TARGETS}: test-default-init-%: KokkosCore_UnitTest_DefaultDeviceTypeInit_%
	./KokkosCore_UnitTest_DefaultDeviceTypeInit_$*

build_all: $(TARGETS)

test: $(TEST_TARGETS)

clean: kokkos-clean
	rm -f *.o $(TARGETS)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(TEST_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

gtest-all.o:$(GTEST_PATH)/gtest/gtest-all.cc
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $(GTEST_PATH)/gtest/gtest-all.cc
