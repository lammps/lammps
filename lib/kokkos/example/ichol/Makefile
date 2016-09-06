SCOTCH_PATH = /home/hcedwar/scotch/6.0.0
KOKKOS_PATH = ../..

vpath %.cpp ${KOKKOS_PATH}/example/ichol/src ${KOKKOS_PATH}/example/ichol/example 

EXAMPLE_HEADERS = $(wildcard $(KOKKOS_PATH)/example/ichol/src/*.hpp ${KOKKOS_PATH}/example/ichol/example/*.hpp )

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

KOKKOS_CXXFLAGS +=	\
	-I${KOKKOS_PATH}/example/ichol/src	\
	-I${KOKKOS_PATH}/example/ichol/example	\
	-I${SCOTCH_PATH}/include

EXE_EXAMPLE_ICHOL_THREADS = KokkosExample_ichol_threads
OBJ_EXAMPLE_ICHOL_THREADS = example_chol_performance_device_pthread.o

EXE_EXAMPLE_ICHOL_CUDA = KokkosExample_ichol_cuda
OBJ_EXAMPLE_ICHOL_CUDA = example_chol_performance_device_cuda.o

TARGETS = $(EXE_EXAMPLE_ICHOL_THREADS) $(EXE_EXAMPLE_ICHOL_CUDA)

#TEST_TARGETS =

$(EXE_EXAMPLE_ICHOL_THREADS) : $(OBJ_EXAMPLE_ICHOL_THREADS) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) \
	$(OBJ_EXAMPLE_ICHOL_THREADS) $(KOKKOS_LIBS) $(LIB) \
	-L${SCOTCH_PATH}/lib -lscotch  -lscotcherr  -lscotcherrexit \
	-o $(EXE_EXAMPLE_ICHOL_THREADS)

$(EXE_EXAMPLE_ICHOL_CUDA) : $(OBJ_EXAMPLE_ICHOL_CUDA) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) \
	$(OBJ_EXAMPLE_ICHOL_CUDA) $(KOKKOS_LIBS) $(LIB) \
	-L${SCOTCH_PATH}/lib -lscotch  -lscotcherr  -lscotcherrexit \
	-o $(EXE_EXAMPLE_ICHOL_CUDA)

build_all : $(TARGETS)

test : build_all

clean: kokkos-clean
	rm -f *.o $(TARGETS)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(EXAMPLE_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

