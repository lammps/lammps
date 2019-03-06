KOKKOS_PATH = ../..

GTEST_PATH = ../../TPL/gtest

vpath %.cpp ${KOKKOS_PATH}/containers/unit_tests
vpath %.cpp ${KOKKOS_PATH}/containers/unit_tests/openmp
vpath %.cpp ${KOKKOS_PATH}/containers/unit_tests/serial
vpath %.cpp ${KOKKOS_PATH}/containers/unit_tests/threads
vpath %.cpp ${KOKKOS_PATH}/containers/unit_tests/rocm
vpath %.cpp ${KOKKOS_PATH}/containers/unit_tests/cuda

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

KOKKOS_CXXFLAGS += -I$(GTEST_PATH) -I${KOKKOS_PATH}/containers/unit_tests

TEST_TARGETS =
TARGETS =

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
	OBJ_CUDA = UnitTestMain.o gtest-all.o
	OBJ_CUDA += TestCuda_BitSet.o
	OBJ_CUDA += TestCuda_DualView.o
	OBJ_CUDA += TestCuda_DynamicView.o
	OBJ_CUDA += TestCuda_DynRankViewAPI_generic.o
	OBJ_CUDA += TestCuda_DynRankViewAPI_rank12345.o
	OBJ_CUDA += TestCuda_DynRankViewAPI_rank67.o
	OBJ_CUDA += TestCuda_ErrorReporter.o
	OBJ_CUDA += TestCuda_OffsetView.o
	OBJ_CUDA += TestCuda_ScatterView.o
	OBJ_CUDA += TestCuda_StaticCrsGraph.o
	OBJ_CUDA += TestCuda_UnorderedMap.o
	OBJ_CUDA += TestCuda_Vector.o
	OBJ_CUDA += TestCuda_ViewCtorPropEmbeddedDim.o
	TARGETS += KokkosContainers_UnitTest_Cuda
	TEST_TARGETS += test-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_ROCM), 1)
	OBJ_ROCM = UnitTestMain.o gtest-all.o
	OBJ_ROCM += TestROCm_BitSet.o
	OBJ_ROCM += TestROCm_DualView.o
	OBJ_ROCM += TestROCm_DynamicView.o
	OBJ_ROCM += TestROCm_DynRankViewAPI_generic.o
	OBJ_ROCM += TestROCm_DynRankViewAPI_rank12345.o
	OBJ_ROCM += TestROCm_DynRankViewAPI_rank67.o
	OBJ_ROCM += TestROCm_ErrorReporter.o
	OBJ_ROCM += TestROCm_OffsetView.o
	OBJ_ROCM += TestROCm_ScatterView.o
	OBJ_ROCM += TestROCm_StaticCrsGraph.o
	OBJ_ROCM += TestROCm_UnorderedMap.o
	OBJ_ROCM += TestROCm_Vector.o
	OBJ_ROCM += TestROCm_ViewCtorPropEmbeddedDim.o
	TARGETS += KokkosContainers_UnitTest_ROCm
	TEST_TARGETS += test-rocm
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)
	OBJ_THREADS = UnitTestMain.o gtest-all.o
	OBJ_THREADS += TestThreads_BitSet.o
	OBJ_THREADS += TestThreads_DualView.o
	OBJ_THREADS += TestThreads_DynamicView.o
	OBJ_THREADS += TestThreads_DynRankViewAPI_generic.o
	OBJ_THREADS += TestThreads_DynRankViewAPI_rank12345.o
	OBJ_THREADS += TestThreads_DynRankViewAPI_rank67.o
	OBJ_THREADS += TestThreads_ErrorReporter.o
	OBJ_THREADS += TestThreads_OffsetView.o
	OBJ_THREADS += TestThreads_ScatterView.o
	OBJ_THREADS += TestThreads_StaticCrsGraph.o
	OBJ_THREADS += TestThreads_UnorderedMap.o
	OBJ_THREADS += TestThreads_Vector.o
	OBJ_THREADS += TestThreads_ViewCtorPropEmbeddedDim.o
	TARGETS += KokkosContainers_UnitTest_Threads
	TEST_TARGETS += test-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
	OBJ_OPENMP = UnitTestMain.o gtest-all.o
	OBJ_OPENMP += TestOpenMP_BitSet.o
	OBJ_OPENMP += TestOpenMP_DualView.o
	OBJ_OPENMP += TestOpenMP_DynamicView.o
	OBJ_OPENMP += TestOpenMP_DynRankViewAPI_generic.o
	OBJ_OPENMP += TestOpenMP_DynRankViewAPI_rank12345.o
	OBJ_OPENMP += TestOpenMP_DynRankViewAPI_rank67.o
	OBJ_OPENMP += TestOpenMP_ErrorReporter.o
	OBJ_OPENMP += TestOpenMP_OffsetView.o
	OBJ_OPENMP += TestOpenMP_ScatterView.o
	OBJ_OPENMP += TestOpenMP_StaticCrsGraph.o
	OBJ_OPENMP += TestOpenMP_UnorderedMap.o
	OBJ_OPENMP += TestOpenMP_Vector.o
	OBJ_OPENMP += TestOpenMP_ViewCtorPropEmbeddedDim.o
	TARGETS += KokkosContainers_UnitTest_OpenMP
	TEST_TARGETS += test-openmp
endif

ifeq ($(KOKKOS_INTERNAL_USE_SERIAL), 1)
	OBJ_SERIAL = UnitTestMain.o gtest-all.o
	OBJ_SERIAL += TestSerial_BitSet.o
	OBJ_SERIAL += TestSerial_DualView.o
	OBJ_SERIAL += TestSerial_DynamicView.o
	OBJ_SERIAL += TestSerial_DynRankViewAPI_generic.o
	OBJ_SERIAL += TestSerial_DynRankViewAPI_rank12345.o
	OBJ_SERIAL += TestSerial_DynRankViewAPI_rank67.o
	OBJ_SERIAL += TestSerial_ErrorReporter.o
	OBJ_SERIAL += TestSerial_OffsetView.o
	OBJ_SERIAL += TestSerial_ScatterView.o
	OBJ_SERIAL += TestSerial_StaticCrsGraph.o
	OBJ_SERIAL += TestSerial_UnorderedMap.o
	OBJ_SERIAL += TestSerial_Vector.o
	OBJ_SERIAL += TestSerial_ViewCtorPropEmbeddedDim.o
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
