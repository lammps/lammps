KOKKOS_PATH = ../..

GTEST_PATH = ../../tpls/gtest

vpath %.cpp ${KOKKOS_PATH}/core/perf_test

default: build_all
	echo "End Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
  CXX = $(KOKKOS_PATH)/bin/nvcc_wrapper
  KOKKOS_CUDA_OPTIONS=enable_lambda
else
  CXX = g++
endif

CXXFLAGS = -O3 
#CXXFLAGS += -DGENERIC_REDUCER
LINK ?= $(CXX)
LDFLAGS ?=
override LDFLAGS += -lpthread

include $(KOKKOS_PATH)/Makefile.kokkos

KOKKOS_CXXFLAGS += -I$(GTEST_PATH) -I${KOKKOS_PATH}/core/perf_test

TEST_TARGETS =
TARGETS =

#

OBJ_PERF = PerfTestMain.o gtest-all.o
OBJ_PERF += PerfTest_ExecSpacePartitioning.o
OBJ_PERF += PerfTestGramSchmidt.o
OBJ_PERF += PerfTestHexGrad.o
OBJ_PERF += PerfTest_CustomReduction.o
OBJ_PERF += PerfTest_ViewCopy_a123.o PerfTest_ViewCopy_b123.o PerfTest_ViewCopy_c123.o PerfTest_ViewCopy_d123.o
OBJ_PERF += PerfTest_ViewCopy_a45.o PerfTest_ViewCopy_b45.o PerfTest_ViewCopy_c45.o PerfTest_ViewCopy_d45.o
OBJ_PERF += PerfTest_ViewCopy_a6.o PerfTest_ViewCopy_b6.o PerfTest_ViewCopy_c6.o PerfTest_ViewCopy_d6.o
OBJ_PERF += PerfTest_ViewCopy_a7.o PerfTest_ViewCopy_b7.o PerfTest_ViewCopy_c7.o PerfTest_ViewCopy_d7.o
OBJ_PERF += PerfTest_ViewCopy_a8.o PerfTest_ViewCopy_b8.o PerfTest_ViewCopy_c8.o PerfTest_ViewCopy_d8.o
OBJ_PERF += PerfTest_ViewAllocate.o
OBJ_PERF += PerfTest_ViewFill_123.o PerfTest_ViewFill_45.o PerfTest_ViewFill_6.o PerfTest_ViewFill_7.o PerfTest_ViewFill_8.o
OBJ_PERF += PerfTest_ViewResize_123.o PerfTest_ViewResize_45.o PerfTest_ViewResize_6.o PerfTest_ViewResize_7.o PerfTest_ViewResize_8.o
TARGETS += KokkosCore_PerformanceTest
TEST_TARGETS += test-performance

#

OBJ_ATOMICS = test_atomic.o 
TARGETS += KokkosCore_PerformanceTest_Atomics
TEST_TARGETS += test-atomic

#

OBJ_MEMPOOL = test_mempool.o 
TARGETS += KokkosCore_PerformanceTest_Mempool
TEST_TARGETS += test-mempool

#

OBJ_TASKDAG = test_taskdag.o 
TARGETS += KokkosCore_PerformanceTest_TaskDAG
TEST_TARGETS += test-taskdag

#

OBJ_ATOMICS_MINMAX = test_atomic_minmax_simple.o
TARGETS += KokkosCore_PerformanceTest_Atomics_MinMax
TEST_TARGETS += test-atomic-minmax

#

KokkosCore_PerformanceTest: $(OBJ_PERF) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_PERF) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_PerformanceTest

KokkosCore_PerformanceTest_Atomics: $(OBJ_ATOMICS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_ATOMICS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_PerformanceTest_Atomics

KokkosCore_PerformanceTest_Mempool: $(OBJ_MEMPOOL) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_MEMPOOL) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_PerformanceTest_Mempool

KokkosCore_PerformanceTest_TaskDAG: $(OBJ_TASKDAG) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_TASKDAG) $(KOKKOS_LIBS) $(LIB) -o KokkosCore_PerformanceTest_TaskDAG

KokkosCore_PerformanceTest_Atomics_MinMax: $(OBJ_ATOMICS_MINMAX) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(EXTRA_PATH) $(OBJ_ATOMICS_MINMAX) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosCore_PerformanceTest_Atomics_MinMax

test-performance: KokkosCore_PerformanceTest
	./KokkosCore_PerformanceTest

test-atomic: KokkosCore_PerformanceTest_Atomics
	./KokkosCore_PerformanceTest_Atomics

test-mempool: KokkosCore_PerformanceTest_Mempool
	./KokkosCore_PerformanceTest_Mempool

test-taskdag: KokkosCore_PerformanceTest_TaskDAG
	./KokkosCore_PerformanceTest_TaskDAG

test-atomic-minmax: KokkosCore_PerformanceTest_Atomics_MinMax
	./KokkosCore_PerformanceTest_Atomics_MinMax

build_all: $(TARGETS)

test: $(TEST_TARGETS)

clean: kokkos-clean
	rm -f *.o $(TARGETS)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

gtest-all.o:$(GTEST_PATH)/gtest/gtest-all.cc
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $(GTEST_PATH)/gtest/gtest-all.cc

