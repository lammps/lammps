KOKKOS_PATH = ../..

GTEST_PATH = ../../TPL/gtest

vpath %.cpp ${KOKKOS_PATH}/algorithms/unit_tests

default: build_all
	echo "End Build"
	

include $(KOKKOS_PATH)/Makefile.kokkos

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
	CXX = nvcc_wrapper
	CXXFLAGS ?= -O3
	LINK = $(CXX)
	LDFLAGS ?= -lpthread
else
	CXX ?= g++
	CXXFLAGS ?= -O3
	LINK ?= $(CXX)
	LDFLAGS ?= -lpthread
endif

KOKKOS_CXXFLAGS += -I$(GTEST_PATH) -I${KOKKOS_PATH}/algorithms/unit_tests

TEST_TARGETS = 
TARGETS = 

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
	OBJ_CUDA = TestCuda.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosAlgorithms_UnitTest_Cuda
	TEST_TARGETS += test-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)
	OBJ_THREADS = TestThreads.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosAlgorithms_UnitTest_Threads
	TEST_TARGETS += test-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
	OBJ_OPENMP = TestOpenMP.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosAlgorithms_UnitTest_OpenMP
	TEST_TARGETS += test-openmp
endif

ifeq ($(KOKKOS_INTERNAL_USE_SERIAL), 1)
	OBJ_SERIAL = TestSerial.o UnitTestMain.o gtest-all.o
	TARGETS += KokkosAlgorithms_UnitTest_Serial
	TEST_TARGETS += test-serial
endif

KokkosAlgorithms_UnitTest_Cuda: $(OBJ_CUDA) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_CUDA) $(KOKKOS_LIBS) $(LIB) -o KokkosAlgorithms_UnitTest_Cuda

KokkosAlgorithms_UnitTest_Threads: $(OBJ_THREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_THREADS) $(KOKKOS_LIBS) $(LIB) -o KokkosAlgorithms_UnitTest_Threads
	
KokkosAlgorithms_UnitTest_OpenMP: $(OBJ_OPENMP) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_OPENMP) $(KOKKOS_LIBS) $(LIB) -o KokkosAlgorithms_UnitTest_OpenMP

KokkosAlgorithms_UnitTest_Serial: $(OBJ_SERIAL) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_SERIAL) $(KOKKOS_LIBS) $(LIB) -o KokkosAlgorithms_UnitTest_Serial

test-cuda: KokkosAlgorithms_UnitTest_Cuda
	./KokkosAlgorithms_UnitTest_Cuda

test-threads: KokkosAlgorithms_UnitTest_Threads
	./KokkosAlgorithms_UnitTest_Threads

test-openmp: KokkosAlgorithms_UnitTest_OpenMP
	./KokkosAlgorithms_UnitTest_OpenMP

test-serial: KokkosAlgorithms_UnitTest_Serial
	./KokkosAlgorithms_UnitTest_Serial
	
build_all: $(TARGETS)

test: $(TEST_TARGETS)
	
clean: kokkos-clean 
	rm -f *.o $(TARGETS)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

gtest-all.o:$(GTEST_PATH)/gtest/gtest-all.cc 
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $(GTEST_PATH)/gtest/gtest-all.cc

