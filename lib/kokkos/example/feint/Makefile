KOKKOS_PATH = ../..

vpath %.cpp ${KOKKOS_PATH}/example/fixture ${KOKKOS_PATH}/example/feint

EXAMPLE_HEADERS = $(wildcard $(KOKKOS_PATH)/example/common/*.hpp ${KOKKOS_PATH}/example/fixture/*.hpp ${KOKKOS_PATH}/example/feint/*.hpp)

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

KOKKOS_CXXFLAGS +=	\
	-I${KOKKOS_PATH}/example/common	\
	-I${KOKKOS_PATH}/example/fixture	\
	-I${KOKKOS_PATH}/example/feint


EXE_EXAMPLE_FEINT = KokkosExample_Feint
OBJ_EXAMPLE_FEINT = BoxElemPart.o main.o

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
  OBJ_EXAMPLE_FEINT += feint_cuda.o
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)
  OBJ_EXAMPLE_FEINT += feint_threads.o
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
  OBJ_EXAMPLE_FEINT += feint_openmp.o
endif

TARGETS = $(EXE_EXAMPLE_FEINT)

#TEST_TARGETS =

$(EXE_EXAMPLE_FEINT) : $(OBJ_EXAMPLE_FEINT) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_EXAMPLE_FEINT) $(KOKKOS_LIBS) $(LIB) -o $(EXE_EXAMPLE_FEINT)

build_all : $(TARGETS)


test : build_all


# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(EXAMPLE_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

