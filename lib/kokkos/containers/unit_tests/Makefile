KOKKOS_PATH = ../..

GTEST_PATH = ../../TPL/gtest

vpath %.cpp ${KOKKOS_PATH}/containers/unit_tests

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

KOKKOS_CXXFLAGS += -I$(GTEST_PATH) -I${KOKKOS_PATH}/containers/unit_tests

TEST_TARGETS =
TARGETS =

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
	OBJ_CUDA = TestCuda.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosContainers_UnitTest_Cuda
	TEST_TARGETS += test-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_ROCM), 1)
	OBJ_ROCM = TestROCm.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosContainers_UnitTest_ROCm
	TEST_TARGETS += test-rocm
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)
	OBJ_THREADS = TestThreads.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosContainers_UnitTest_Threads
	TEST_TARGETS += test-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
	OBJ_OPENMP = TestOpenMP.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosContainers_UnitTest_OpenMP
	TEST_TARGETS += test-openmp
endif

ifeq ($(KOKKOS_INTERNAL_USE_SERIAL), 1)
	OBJ_SERIAL = TestSerial.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosContainers_UnitTest_Serial
	TEST_TARGETS += test-serial
endif

KokkosContainers_UnitTest_Cuda: $(OBJ_CUDA) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_CUDA) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosContainers_UnitTest_Cuda

KokkosContainers_UnitTest_ROCm: $(OBJ_ROCM) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_ROCM) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosContainers_UnitTest_ROCm

KokkosContainers_UnitTest_Threads: $(OBJ_THREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_THREADS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosContainers_UnitTest_Threads

KokkosContainers_UnitTest_OpenMP: $(OBJ_OPENMP) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_OPENMP) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosContainers_UnitTest_OpenMP

KokkosContainers_UnitTest_Serial: $(OBJ_SERIAL) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_SERIAL) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosContainers_UnitTest_Serial

test-cuda: KokkosContainers_UnitTest_Cuda
	./KokkosContainers_UnitTest_Cuda

test-rocm: KokkosContainers_UnitTest_ROCm
	./KokkosContainers_UnitTest_ROCm

test-threads: KokkosContainers_UnitTest_Threads
	./KokkosContainers_UnitTest_Threads

test-openmp: KokkosContainers_UnitTest_OpenMP
	./KokkosContainers_UnitTest_OpenMP

test-serial: KokkosContainers_UnitTest_Serial
	./KokkosContainers_UnitTest_Serial

build_all: $(TARGETS)

test: $(TEST_TARGETS)

clean: kokkos-clean
	rm -f *.o $(TARGETS)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

gtest-all.o:$(GTEST_PATH)/gtest/gtest-all.cc
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $(GTEST_PATH)/gtest/gtest-all.cc
