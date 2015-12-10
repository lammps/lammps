KOKKOS_PATH = ../..

GTEST_PATH = ../../tpls/gtest

vpath %.cpp ${KOKKOS_PATH}/core/unit_test
TEST_HEADERS = $(wildcard $(KOKKOS_PATH)/core/unit_test/*.hpp)

default: build_all
	echo "End Build"

include $(KOKKOS_PATH)/Makefile.kokkos

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
	CXX = $(NVCC_WRAPPER)
	CXXFLAGS ?= -O3
	LINK = $(CXX)
	LDFLAGS ?= -lpthread
else
	CXX ?= g++
	CXXFLAGS ?= -O3
	LINK ?= $(CXX)
	LDFLAGS ?= -lpthread
endif

KOKKOS_CXXFLAGS += -I$(GTEST_PATH) -I${KOKKOS_PATH}/core/unit_test

TEST_TARGETS =
TARGETS =

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
	OBJ_CUDA = TestCuda.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_Cuda
	TEST_TARGETS += test-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)
	OBJ_THREADS = TestThreads.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_Threads
	TEST_TARGETS += test-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
	OBJ_OPENMP = TestOpenMP.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_OpenMP
	TEST_TARGETS += test-openmp
endif

ifeq ($(KOKKOS_INTERNAL_USE_SERIAL), 1)
	OBJ_SERIAL = TestSerial.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_Serial
	TEST_TARGETS += test-serial
endif

ifeq ($(KOKKOS_INTERNAL_USE_QTHREAD), 1)
	OBJ_QTHREAD = TestQthread.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosCore_UnitTest_Qthread
	TEST_TARGETS += test-qthread
endif

OBJ_HWLOC = TestHWLOC.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_HWLOC
TEST_TARGETS += test-hwloc

OBJ_ALLOCATIONTRACKER = TestAllocationTracker.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_AllocationTracker
TEST_TARGETS += test-allocationtracker

OBJ_DEFAULT = TestDefaultDeviceType.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_Default
TEST_TARGETS += test-default

OBJ_DEFAULTINIT = TestDefaultDeviceTypeInit.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_DefaultInit
TEST_TARGETS += test-default-init

OBJ_SYNCHRONIC = TestSynchronic.o UnitTestMain.o gtest-all.o
TARGETS += KokkosCore_UnitTest_Synchronic
TEST_TARGETS += test-synchronic

KokkosCore_UnitTest_Cuda: $(OBJ_CUDA) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_CUDA) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Cuda

KokkosCore_UnitTest_Threads: $(OBJ_THREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_THREADS) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Threads

KokkosCore_UnitTest_OpenMP: $(OBJ_OPENMP) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_OPENMP) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_OpenMP

KokkosCore_UnitTest_Serial: $(OBJ_SERIAL) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_SERIAL) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Serial

KokkosCore_UnitTest_Qthread: $(OBJ_QTHREAD) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_QTHREAD) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Qthread

KokkosCore_UnitTest_HWLOC: $(OBJ_HWLOC) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_HWLOC) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_HWLOC

KokkosCore_UnitTest_AllocationTracker: $(OBJ_ALLOCATIONTRACKER) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_ALLOCATIONTRACKER) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_AllocationTracker

KokkosCore_UnitTest_Default: $(OBJ_DEFAULT) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_DEFAULT) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Default

KokkosCore_UnitTest_DefaultInit: $(OBJ_DEFAULTINIT) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_DEFAULTINIT) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_DefaultInit

KokkosCore_UnitTest_Synchronic: $(OBJ_SYNCHRONIC) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_SYNCHRONIC) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_UnitTest_Synchronic

test-cuda: KokkosCore_UnitTest_Cuda
	./KokkosCore_UnitTest_Cuda

test-threads: KokkosCore_UnitTest_Threads
	./KokkosCore_UnitTest_Threads

test-openmp: KokkosCore_UnitTest_OpenMP
	./KokkosCore_UnitTest_OpenMP

test-serial: KokkosCore_UnitTest_Serial
	./KokkosCore_UnitTest_Serial

test-qthread: KokkosCore_UnitTest_Qthread
	./KokkosCore_UnitTest_Qthread

test-hwloc: KokkosCore_UnitTest_HWLOC
	./KokkosCore_UnitTest_HWLOC

test-allocationtracker: KokkosCore_UnitTest_AllocationTracker
	./KokkosCore_UnitTest_AllocationTracker

test-default: KokkosCore_UnitTest_Default
	./KokkosCore_UnitTest_Default

test-default-init: KokkosCore_UnitTest_DefaultInit
	./KokkosCore_UnitTest_DefaultInit

test-synchronic: KokkosCore_UnitTest_Synchronic
	./KokkosCore_UnitTest_Synchronic

build_all: $(TARGETS)

test: $(TEST_TARGETS)

clean: kokkos-clean
	rm -f *.o $(TARGETS)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(TEST_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

gtest-all.o:$(GTEST_PATH)/gtest/gtest-all.cc
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $(GTEST_PATH)/gtest/gtest-all.cc

