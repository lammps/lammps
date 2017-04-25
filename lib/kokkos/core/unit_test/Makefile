KOKKOS_PATH = ../..

GTEST_PATH = ../../tpls/gtest

vpath %.cpp ${KOKKOS_PATH}/core/unit_test
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/serial
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/threads
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/openmp
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/qthreads
vpath %.cpp ${KOKKOS_PATH}/core/unit_test/cuda

TEST_HEADERS = $(wildcard $(KOKKOS_PATH)/core/unit_test/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOS_PATH)/core/unit_test/*/*.hpp)

default: build_all
	echo "End Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
  CXX = $(KOKKOS_PATH)/config/nvcc_wrapper
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
	OBJ_CUDA = TestCuda_Other.o TestCuda_Reductions_a.o TestCuda_Reductions_b.o TestCuda_Atomics.o TestCuda_Team.o TestCuda_Spaces.o
	OBJ_CUDA += TestCuda_SubView_a.o TestCuda_SubView_b.o
ifeq ($(KOKKOS_INTERNAL_COMPILER_XL), 1)
	OBJ_OPENMP += TestCuda_SubView_c_all.o
else
	OBJ_CUDA += TestCuda_SubView_c01.o TestCuda_SubView_c02.o TestCuda_SubView_c03.o
	OBJ_CUDA += TestCuda_SubView_c04.o TestCuda_SubView_c05.o TestCuda_SubView_c06.o
	OBJ_CUDA += TestCuda_SubView_c07.o TestCuda_SubView_c08.o TestCuda_SubView_c09.o
	OBJ_CUDA += TestCuda_SubView_c10.o TestCuda_SubView_c11.o TestCuda_SubView_c12.o
endif
	OBJ_CUDA += TestCuda_ViewAPI_a.o TestCuda_ViewAPI_b.o TestCuda_ViewAPI_c.o TestCuda_ViewAPI_d.o
	OBJ_CUDA += TestCuda_ViewAPI_e.o TestCuda_ViewAPI_f.o TestCuda_ViewAPI_g.o TestCuda_ViewAPI_h.o
	OBJ_CUDA += TestCuda_ViewAPI_s.o
	OBJ_CUDA += UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_Cuda
	TEST_TARGETS += test-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)
	OBJ_THREADS = TestThreads_Other.o TestThreads_Reductions.o TestThreads_Atomics.o TestThreads_Team.o
	OBJ_THREADS += TestThreads_SubView_a.o TestThreads_SubView_b.o
	OBJ_THREADS += TestThreads_SubView_c01.o TestThreads_SubView_c02.o TestThreads_SubView_c03.o
	OBJ_THREADS += TestThreads_SubView_c04.o TestThreads_SubView_c05.o TestThreads_SubView_c06.o
	OBJ_THREADS += TestThreads_SubView_c07.o TestThreads_SubView_c08.o TestThreads_SubView_c09.o
	OBJ_THREADS += TestThreads_SubView_c10.o TestThreads_SubView_c11.o TestThreads_SubView_c12.o
	OBJ_THREADS += TestThreads_ViewAPI_a.o TestThreads_ViewAPI_b.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_Threads
	TEST_TARGETS += test-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
	OBJ_OPENMP = TestOpenMP_Other.o TestOpenMP_Reductions.o TestOpenMP_Atomics.o TestOpenMP_Team.o
	OBJ_OPENMP += TestOpenMP_SubView_a.o TestOpenMP_SubView_b.o
ifeq ($(KOKKOS_INTERNAL_COMPILER_XL), 1)
	OBJ_OPENMP += TestOpenMP_SubView_c_all.o
else
	OBJ_OPENMP += TestOpenMP_SubView_c01.o TestOpenMP_SubView_c02.o TestOpenMP_SubView_c03.o
	OBJ_OPENMP += TestOpenMP_SubView_c04.o TestOpenMP_SubView_c05.o TestOpenMP_SubView_c06.o
	OBJ_OPENMP += TestOpenMP_SubView_c07.o TestOpenMP_SubView_c08.o TestOpenMP_SubView_c09.o
	OBJ_OPENMP += TestOpenMP_SubView_c10.o TestOpenMP_SubView_c11.o TestOpenMP_SubView_c12.o
endif
	OBJ_OPENMP += TestOpenMP_ViewAPI_a.o TestOpenMP_ViewAPI_b.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_OpenMP
	TEST_TARGETS += test-openmp
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
	TEST_TARGETS += test-qthreads
endif

ifeq ($(KOKKOS_INTERNAL_USE_SERIAL), 1)
	OBJ_SERIAL = TestSerial_Other.o TestSerial_Reductions.o TestSerial_Atomics.o TestSerial_Team.o
	OBJ_SERIAL += TestSerial_SubView_a.o TestSerial_SubView_b.o
ifeq ($(KOKKOS_INTERNAL_COMPILER_XL), 1)
	OBJ_OPENMP += TestSerial_SubView_c_all.o
else
	OBJ_SERIAL += TestSerial_SubView_c01.o TestSerial_SubView_c02.o TestSerial_SubView_c03.o
	OBJ_SERIAL += TestSerial_SubView_c04.o TestSerial_SubView_c05.o TestSerial_SubView_c06.o
	OBJ_SERIAL += TestSerial_SubView_c07.o TestSerial_SubView_c08.o TestSerial_SubView_c09.o
	OBJ_SERIAL += TestSerial_SubView_c10.o TestSerial_SubView_c11.o TestSerial_SubView_c12.o
endif
	OBJ_SERIAL += TestSerial_ViewAPI_a.o TestSerial_ViewAPI_b.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_Serial
	TEST_TARGETS += test-serial
endif

OBJ_HWLOC = TestHWLOC.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_HWLOC
TEST_TARGETS += test-hwloc

OBJ_DEFAULT = TestDefaultDeviceType.o TestDefaultDeviceType_a.o TestDefaultDeviceType_b.o TestDefaultDeviceType_c.o TestDefaultDeviceType_d.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_Default
TEST_TARGETS += test-default

NUM_INITTESTS = 16
INITTESTS_NUMBERS := $(shell seq 1 ${NUM_INITTESTS})
INITTESTS_TARGETS := $(addprefix KokkosCore_UnitTest_DefaultDeviceTypeInit_,${INITTESTS_NUMBERS})
TARGETS += ${INITTESTS_TARGETS}
INITTESTS_TEST_TARGETS := $(addprefix test-default-init-,${INITTESTS_NUMBERS})
TEST_TARGETS += ${INITTESTS_TEST_TARGETS}

KokkosCore_UnitTest_Cuda: $(OBJ_CUDA) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_CUDA) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Cuda

KokkosCore_UnitTest_Threads: $(OBJ_THREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_THREADS) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Threads

KokkosCore_UnitTest_OpenMP: $(OBJ_OPENMP) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_OPENMP) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_OpenMP

KokkosCore_UnitTest_Serial: $(OBJ_SERIAL) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_SERIAL) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Serial

KokkosCore_UnitTest_Qthreads: $(OBJ_QTHREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_QTHREADS) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Qthreads

KokkosCore_UnitTest_HWLOC: $(OBJ_HWLOC) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_HWLOC) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_HWLOC

KokkosCore_UnitTest_AllocationTracker: $(OBJ_ALLOCATIONTRACKER) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_ALLOCATIONTRACKER) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_AllocationTracker

KokkosCore_UnitTest_Default: $(OBJ_DEFAULT) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_DEFAULT) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Default

${INITTESTS_TARGETS}: KokkosCore_UnitTest_DefaultDeviceTypeInit_%: TestDefaultDeviceTypeInit_%.o UnitTestMain.o gtest-all.o $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) TestDefaultDeviceTypeInit_$*.o UnitTestMain.o gtest-all.o $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_DefaultDeviceTypeInit_$*

test-cuda: KokkosCore_UnitTest_Cuda
	./KokkosCore_UnitTest_Cuda

test-threads: KokkosCore_UnitTest_Threads
	./KokkosCore_UnitTest_Threads

test-openmp: KokkosCore_UnitTest_OpenMP
	./KokkosCore_UnitTest_OpenMP

test-serial: KokkosCore_UnitTest_Serial
	./KokkosCore_UnitTest_Serial

test-qthreads: KokkosCore_UnitTest_Qthreads
	./KokkosCore_UnitTest_Qthreads

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
