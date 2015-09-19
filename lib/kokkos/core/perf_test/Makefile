KOKKOS_PATH = ../..

GTEST_PATH = ../../TPL/gtest

vpath %.cpp ${KOKKOS_PATH}/core/perf_test

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
	LDFLAGS ?=  -lpthread
endif

KOKKOS_CXXFLAGS += -I$(GTEST_PATH) -I${KOKKOS_PATH}/core/perf_test

TEST_TARGETS = 
TARGETS = 

OBJ_PERF = PerfTestHost.o PerfTestCuda.o PerfTestMain.o gtest-all.o
TARGETS += KokkosCore_PerformanceTest
TEST_TARGETS += test-performance

OBJ_ATOMICS = test_atomic.o 
TARGETS += KokkosCore_PerformanceTest_Atomics
TEST_TARGETS += test-atomic


KokkosCore_PerformanceTest: $(OBJ_PERF) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_PERF) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_PerformanceTest

KokkosCore_PerformanceTest_Atomics: $(OBJ_ATOMICS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_ATOMICS) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_PerformanceTest_Atomics

test-performance: KokkosCore_PerformanceTest
	./KokkosCore_PerformanceTest

test-atomic: KokkosCore_PerformanceTest_Atomics
	./KokkosCore_PerformanceTest_Atomics
	

build_all: $(TARGETS)

test: $(TEST_TARGETS)
	
clean: kokkos-clean 
	rm -f *.o $(TARGETS)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

gtest-all.o:$(GTEST_PATH)/gtest/gtest-all.cc 
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $(GTEST_PATH)/gtest/gtest-all.cc

