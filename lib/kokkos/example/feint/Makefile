KOKKOS_PATH = ../..
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
vpath %.cpp ${KOKKOS_SRC_PATH}/example/fixture ${KOKKOS_SRC_PATH}/example/feint

EXAMPLE_HEADERS = $(wildcard $(KOKKOS_SRC_PATH)/example/common/*.hpp ${KOKKOS_SRC_PATH}/example/fixture/*.hpp ${KOKKOS_SRC_PATH}/example/feint/*.hpp)

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

include $(KOKKOS_PATH)/Makefile.kokkos        

KOKKOS_CXXFLAGS +=	\
	-I${KOKKOS_SRC_PATH}/example/common	\
	-I${KOKKOS_SRC_PATH}/example/fixture	\
	-I${KOKKOS_SRC_PATH}/example/feint

EXE_EXAMPLE_FEINT = KokkosExample_Feint
OBJ_EXAMPLE_FEINT = BoxElemPart.o main.o

ifeq ($(KOKKOS_INTERNAL_USE_SERIAL), 1)
  OBJ_EXAMPLE_FEINT += feint_serial.o
endif

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
  OBJ_EXAMPLE_FEINT += feint_cuda.o
endif

ifeq ($(KOKKOS_INTERNAL_USE_ROCM), 1)
  OBJ_EXAMPLE_FEINT += feint_rocm.o
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

clean: kokkos-clean
	rm -f *.o $(TARGETS)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(EXAMPLE_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<
