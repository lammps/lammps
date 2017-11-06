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
LDFLAGS ?= -lpthread

include $(KOKKOS_PATH)/Makefile.kokkos

KOKKOS_CXXFLAGS += -I$(GTEST_PATH) -I${KOKKOS_PATH}/core/unit_test

TEST_TARGETS =
TARGETS =

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
	OBJ_CUDA = UnitTestMainInit.o gtest-all.o
	OBJ_CUDA += TestCuda_Init.o
	OBJ_CUDA += TestCuda_SharedAlloc.o TestCudaUVM_SharedAlloc.o TestCudaHostPinned_SharedAlloc.o
	OBJ_CUDA += TestCuda_RangePolicy.o
	OBJ_CUDA += TestCuda_ViewAPI_b.o
	OBJ_CUDA += TestCuda_ViewMapping_a.o TestCuda_ViewMapping_b.o TestCuda_ViewMapping_subview.o
	OBJ_CUDA += TestCudaUVM_ViewAPI.o
	OBJ_CUDA += TestCudaUVM_ViewMapping_a.o TestCudaUVM_ViewMapping_b.o TestCudaUVM_ViewMapping_subview.o
	OBJ_CUDA += TestCudaHostPinned_ViewAPI.o
	OBJ_CUDA += TestCudaHostPinned_ViewMapping_a.o TestCudaHostPinned_ViewMapping_b.o TestCudaHostPinned_ViewMapping_subview.o
	OBJ_CUDA += TestCuda_ViewOfClass.o
ifeq ($(KOKKOS_INTERNAL_COMPILER_XL), 1)
	OBJ_CUDA += TestCuda_SubView_c_all.o
else
	OBJ_CUDA += TestCuda_SubView_a.o TestCuda_SubView_b.o
	OBJ_CUDA += TestCuda_SubView_c01.o TestCuda_SubView_c02.o TestCuda_SubView_c03.o
	OBJ_CUDA += TestCuda_SubView_c04.o TestCuda_SubView_c05.o TestCuda_SubView_c06.o
	OBJ_CUDA += TestCuda_SubView_c07.o TestCuda_SubView_c08.o TestCuda_SubView_c09.o
	OBJ_CUDA += TestCuda_SubView_c10.o TestCuda_SubView_c11.o TestCuda_SubView_c12.o
	OBJ_CUDA += TestCuda_SubView_c13.o
endif
	OBJ_CUDA += TestCuda_Reductions.o TestCuda_Scan.o
	OBJ_CUDA += TestCuda_Complex.o
	OBJ_CUDA += TestCuda_AtomicOperations.o TestCuda_AtomicViews.o TestCuda_Atomics.o
	OBJ_CUDA += TestCuda_Team.o TestCuda_TeamScratch.o
	OBJ_CUDA += TestCuda_TeamReductionScan.o
	OBJ_CUDA += TestCuda_Other.o
	OBJ_CUDA += TestCuda_MDRange.o
	OBJ_CUDA += TestCuda_Crs.o
	OBJ_CUDA += TestCuda_Task.o TestCuda_WorkGraph.o
	OBJ_CUDA += TestCuda_Spaces.o
	OBJ_CUDA += TestCuda_UniqueToken.o
	
	TARGETS += KokkosCore_UnitTest_Cuda

	TEST_TARGETS += test-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_ROCM), 1)
        OBJ_ROCM = UnitTestMainInit.o gtest-all.o
        OBJ_ROCM += TestROCm_Init.o
        OBJ_ROCM += TestROCm_Complex.o
#        OBJ_ROCM += TestROCm_RangePolicy.o
# rocm.range_scan locking up
        OBJ_ROCM += TestROCm_AtomicOperations.o
        OBJ_ROCM += TestROCm_Atomics.o
# complex failing
        OBJ_ROCM += TestROCm_AtomicViews.o
        OBJ_ROCM += TestROCm_Other.o
# rocm.memory_pool
        OBJ_ROCM += TestROCm_Scan.o
        OBJ_ROCM += TestROCm_SharedAlloc.o
        OBJ_ROCM += TestROCm_SubView_a.o
#        OBJ_ROCM += TestROCm_SubView_b.o
# relies on host accessable device memory
#        OBJ_ROCM += TestROCm_SubView_c01.o
#        OBJ_ROCM += TestROCm_SubView_c02.o
#        OBJ_ROCM += TestROCm_SubView_c03.o
#        OBJ_ROCM += TestROCm_SubView_c04.o
#        OBJ_ROCM += TestROCm_SubView_c05.o
#        OBJ_ROCM += TestROCm_SubView_c06.o
#        OBJ_ROCM += TestROCm_SubView_c07.o
#        OBJ_ROCM += TestROCm_SubView_c08.o
#        OBJ_ROCM += TestROCm_SubView_c09.o
#        OBJ_ROCM += TestROCm_SubView_c10.o
#        OBJ_ROCM += TestROCm_SubView_c11.o
#        OBJ_ROCM += TestROCm_SubView_c12.o
# all of the above use UVM or Host accessable memory
#        OBJ_ROCM += TestROCm_Team.o
# compile fails
#        OBJ_ROCM += TestROCm_TeamReductionScan.o
# compile fails
#        OBJ_ROCM += TestROCm_TeamScratch.o
# compile fails
        OBJ_ROCM += TestROCm_ViewAPI_b.o
#  test fail in view_api
        OBJ_ROCM += TestROCm_ViewMapping_a.o
        OBJ_ROCM += TestROCm_ViewMapping_b.o
        OBJ_ROCM += TestROCm_ViewMapping_subview.o
	OBJ_ROCM += TestROCmHostPinned_ViewAPI.o
	OBJ_ROCM += TestROCmHostPinned_ViewMapping_a.o 
	OBJ_ROCM += TestROCmHostPinned_ViewMapping_b.o 
	OBJ_ROCM += TestROCmHostPinned_ViewMapping_subview.o
        OBJ_ROCM += TestROCm_ViewOfClass.o
	OBJ_ROCM += TestROCm_Spaces.o
     
        TARGETS += KokkosCore_UnitTest_ROCm
        TEST_TARGETS += test-rocm
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)

	OBJ_THREADS = UnitTestMainInit.o gtest-all.o
	OBJ_THREADS += TestThreads_Init.o
	OBJ_THREADS += TestThreads_SharedAlloc.o
	OBJ_THREADS += TestThreads_RangePolicy.o
	OBJ_THREADS += TestThreads_ViewAPI_b.o
	OBJ_THREADS += TestThreads_ViewMapping_a.o TestThreads_ViewMapping_b.o TestThreads_ViewMapping_subview.o
	OBJ_THREADS += TestThreads_ViewOfClass.o
	OBJ_THREADS += TestThreads_SubView_a.o TestThreads_SubView_b.o
	OBJ_THREADS += TestThreads_SubView_c01.o TestThreads_SubView_c02.o TestThreads_SubView_c03.o
	OBJ_THREADS += TestThreads_SubView_c04.o TestThreads_SubView_c05.o TestThreads_SubView_c06.o
	OBJ_THREADS += TestThreads_SubView_c07.o TestThreads_SubView_c08.o TestThreads_SubView_c09.o
	OBJ_THREADS += TestThreads_SubView_c10.o TestThreads_SubView_c11.o TestThreads_SubView_c12.o
	OBJ_THREADS += TestThreads_Reductions.o TestThreads_Scan.o
	OBJ_THREADS += TestThreads_Complex.o
	OBJ_THREADS += TestThreads_AtomicOperations.o TestThreads_AtomicViews.o TestThreads_Atomics.o
	OBJ_THREADS += TestThreads_Team.o TestThreads_TeamScratch.o
	OBJ_THREADS += TestThreads_TeamReductionScan.o
	OBJ_THREADS += TestThreads_Other.o
	OBJ_THREADS += TestThreads_MDRange.o

	TARGETS += KokkosCore_UnitTest_Threads

	TEST_TARGETS += test-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
	OBJ_OPENMP = UnitTestMainInit.o gtest-all.o
	OBJ_OPENMP += TestOpenMP_Init.o
	OBJ_OPENMP += TestOpenMP_SharedAlloc.o
	OBJ_OPENMP += TestOpenMP_RangePolicy.o
	OBJ_OPENMP += TestOpenMP_ViewAPI_b.o
	OBJ_OPENMP += TestOpenMP_ViewMapping_a.o TestOpenMP_ViewMapping_b.o TestOpenMP_ViewMapping_subview.o
	OBJ_OPENMP += TestOpenMP_ViewOfClass.o
ifeq ($(KOKKOS_INTERNAL_COMPILER_XL), 1)
	OBJ_OPENMP += TestOpenMP_SubView_c_all.o
else
	OBJ_OPENMP += TestOpenMP_SubView_a.o TestOpenMP_SubView_b.o
	OBJ_OPENMP += TestOpenMP_SubView_c01.o TestOpenMP_SubView_c02.o TestOpenMP_SubView_c03.o
	OBJ_OPENMP += TestOpenMP_SubView_c04.o TestOpenMP_SubView_c05.o TestOpenMP_SubView_c06.o
	OBJ_OPENMP += TestOpenMP_SubView_c07.o TestOpenMP_SubView_c08.o TestOpenMP_SubView_c09.o
	OBJ_OPENMP += TestOpenMP_SubView_c10.o TestOpenMP_SubView_c11.o TestOpenMP_SubView_c12.o
	OBJ_OPENMP += TestOpenMP_SubView_c13.o
endif
	OBJ_OPENMP += TestOpenMP_Reductions.o TestOpenMP_Scan.o
	OBJ_OPENMP += TestOpenMP_Complex.o
	OBJ_OPENMP += TestOpenMP_AtomicOperations.o TestOpenMP_AtomicViews.o TestOpenMP_Atomics.o
	OBJ_OPENMP += TestOpenMP_Team.o TestOpenMP_TeamScratch.o
	OBJ_OPENMP += TestOpenMP_TeamReductionScan.o
	OBJ_OPENMP += TestOpenMP_Other.o
	OBJ_OPENMP += TestOpenMP_MDRange.o
	OBJ_OPENMP += TestOpenMP_Crs.o
	OBJ_OPENMP += TestOpenMP_Task.o TestOpenMP_WorkGraph.o
	OBJ_OPENMP += TestOpenMP_UniqueToken.o
	
	TARGETS += KokkosCore_UnitTest_OpenMP

	TEST_TARGETS += test-openmp
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMPTARGET), 1)
	OBJ_OPENMPTARGET = UnitTestMainInit.o gtest-all.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_Init.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_SharedAlloc.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_RangePolicy.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_ViewAPI_b.o #Some commented out code
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
        #OBJ_OPENMPTARGET += TestOpenMPTarget_Scan.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_Complex.o
	OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicOperations.o 
        OBJ_OPENMPTARGET += TestOpenMPTarget_AtomicViews.o 
        OBJ_OPENMPTARGET += TestOpenMPTarget_Atomics.o # Commented Out Arbitrary Type Atomics
	#OBJ_OPENMPTARGET += TestOpenMPTarget_Team.o # There is still a static function in this
        #OBJ_OPENMPTARGET += TestOpenMPTarget_TeamScratch.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_TeamReductionScan.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_Other.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_MDRange.o
	#OBJ_OPENMPTARGET += TestOpenMPTarget_Task.o
	
	TARGETS += KokkosCore_UnitTest_OpenMPTarget

	TEST_TARGETS += test-openmptarget

endif

ifeq ($(KOKKOS_INTERNAL_USE_QTHREADS), 1)
	OBJ_QTHREADS = TestQthreads_Other.o TestQthreads_Reductions.o TestQthreads_Atomics.o TestQthreads_Team.o
	OBJ_QTHREADS += TestQthreads_SubView_a.o TestQthreads_SubView_b.o
ifeq ($(KOKKOS_INTERNAL_COMPILER_XL), 1)
	OBJ_QTHREADS += TestQthreads_SubView_c_all.o
else
	OBJ_QTHREADS += TestQthreads_SubView_c01.o TestQthreads_SubView_c02.o TestQthreads_SubView_c03.o
	OBJ_QTHREADS += TestQthreads_SubView_c04.o TestQthreads_SubView_c05.o TestQthreads_SubView_c06.o
	OBJ_QTHREADS += TestQthreads_SubView_c07.o TestQthreads_SubView_c08.o TestQthreads_SubView_c09.o
	OBJ_QTHREADS += TestQthreads_SubView_c10.o TestQthreads_SubView_c11.o TestQthreads_SubView_c12.o
endif
	OBJ_QTHREADS += TestQthreads_ViewAPI_a.o TestQthreads_ViewAPI_b.o UnitTestMain.o gtest-all.o
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
        OBJ_SERIAL += TestSerial_ViewAPI_b.o
        OBJ_SERIAL += TestSerial_ViewMapping_a.o TestSerial_ViewMapping_b.o TestSerial_ViewMapping_subview.o
        OBJ_SERIAL += TestSerial_ViewOfClass.o
ifeq ($(KOKKOS_INTERNAL_COMPILER_XL), 1)
	    OBJ_SERIAL += TestSerial_SubView_c_all.o
else
        OBJ_SERIAL += TestSerial_SubView_a.o TestSerial_SubView_b.o
        OBJ_SERIAL += TestSerial_SubView_c01.o TestSerial_SubView_c02.o TestSerial_SubView_c03.o
        OBJ_SERIAL += TestSerial_SubView_c04.o TestSerial_SubView_c05.o TestSerial_SubView_c06.o
        OBJ_SERIAL += TestSerial_SubView_c07.o TestSerial_SubView_c08.o TestSerial_SubView_c09.o
        OBJ_SERIAL += TestSerial_SubView_c10.o TestSerial_SubView_c11.o TestSerial_SubView_c12.o
        OBJ_SERIAL += TestSerial_SubView_c13.o
endif
        OBJ_SERIAL += TestSerial_Reductions.o TestSerial_Scan.o
        OBJ_SERIAL += TestSerial_Complex.o
        OBJ_SERIAL += TestSerial_AtomicOperations.o TestSerial_AtomicViews.o TestSerial_Atomics.o
        OBJ_SERIAL += TestSerial_Team.o TestSerial_TeamScratch.o
        OBJ_SERIAL += TestSerial_TeamReductionScan.o
        OBJ_SERIAL += TestSerial_Other.o
        #HCC_WORKAROUND
        ifneq ($(KOKKOS_INTERNAL_COMPILER_HCC), 1)
        OBJ_SERIAL += TestSerial_MDRange.o
        endif
        OBJ_SERIAL += TestSerial_Crs.o
        OBJ_SERIAL += TestSerial_Task.o TestSerial_WorkGraph.o
	
	TARGETS += KokkosCore_UnitTest_Serial

	TEST_TARGETS += test-serial
endif

OBJ_HWLOC = TestHWLOC.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_HWLOC
TEST_TARGETS += test-hwloc

OBJ_DEFAULT = UnitTestMainInit.o gtest-all.o
ifneq ($(KOKKOS_INTERNAL_USE_OPENMPTARGET), 1)
ifneq ($(KOKKOS_INTERNAL_COMPILER_HCC), 1)
  OBJ_DEFAULT += TestDefaultDeviceType.o TestDefaultDeviceType_a.o TestDefaultDeviceType_b.o TestDefaultDeviceType_c.o TestDefaultDeviceType_d.o
endif
endif

TARGETS += KokkosCore_UnitTest_Default
TEST_TARGETS += test-default

NUM_INITTESTS = 16
INITTESTS_NUMBERS := $(shell seq 1 ${NUM_INITTESTS})
INITTESTS_TARGETS := $(addprefix KokkosCore_UnitTest_DefaultDeviceTypeInit_,${INITTESTS_NUMBERS})
TARGETS += ${INITTESTS_TARGETS}
INITTESTS_TEST_TARGETS := $(addprefix test-default-init-,${INITTESTS_NUMBERS})
TEST_TARGETS += ${INITTESTS_TEST_TARGETS}

KokkosCore_UnitTest_Cuda: $(OBJ_CUDA) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_CUDA) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Cuda

KokkosCore_UnitTest_ROCm: $(OBJ_ROCM) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_ROCM) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_ROCm

KokkosCore_UnitTest_Threads: $(OBJ_THREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_THREADS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Threads

KokkosCore_UnitTest_OpenMP: $(OBJ_OPENMP) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_OPENMP) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_OpenMP

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

KokkosCore_UnitTest_AllocationTracker: $(OBJ_ALLOCATIONTRACKER) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_ALLOCATIONTRACKER) $(KOKKOS_LIBS) $( $(KOKKOS_LDFLAGS) $(LDFLAGS)LIB) -o KokkosCore_UnitTest_AllocationTracker

KokkosCore_UnitTest_Default: $(OBJ_DEFAULT) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_DEFAULT) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_Default

${INITTESTS_TARGETS}: KokkosCore_UnitTest_DefaultDeviceTypeInit_%: TestDefaultDeviceTypeInit_%.o UnitTestMain.o gtest-all.o $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) TestDefaultDeviceTypeInit_$*.o UnitTestMain.o gtest-all.o $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_UnitTest_DefaultDeviceTypeInit_$*

test-cuda: KokkosCore_UnitTest_Cuda
	./KokkosCore_UnitTest_Cuda

test-rocm: KokkosCore_UnitTest_ROCm
	./KokkosCore_UnitTest_ROCm

test-threads: KokkosCore_UnitTest_Threads
	./KokkosCore_UnitTest_Threads

test-openmp: KokkosCore_UnitTest_OpenMP
	./KokkosCore_UnitTest_OpenMP

test-openmptarget: KokkosCore_UnitTest_OpenMPTarget
	./KokkosCore_UnitTest_OpenMPTarget

test-serial: KokkosCore_UnitTest_Serial
	./KokkosCore_UnitTest_Serial

test-qthreads: KokkosCore_UnitTest_Qthreads KokkosCore_UnitTest_Qthreads2
	./KokkosCore_UnitTest_Qthreads
	./KokkosCore_UnitTest_Qthreads2

test-hwloc: KokkosCore_UnitTest_HWLOC
	./KokkosCore_UnitTest_HWLOC

test-allocationtracker: KokkosCore_UnitTest_AllocationTracker
	./KokkosCore_UnitTest_AllocationTracker

test-default: KokkosCore_UnitTest_Default
	./KokkosCore_UnitTest_Default

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
