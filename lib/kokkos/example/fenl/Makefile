KOKKOS_PATH ?= ../..

MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
SRC_DIR := $(dir $(MAKEFILE_PATH))

vpath %.cpp ${SRC_DIR}/../fixture ${SRC_DIR}

EXAMPLE_HEADERS = $(wildcard $(SRC_DIR)/../common/*.hpp ${SRC_DIR}/../fixture/*.hpp ${SRC_DIR}/*.hpp)

default: build_all
	echo "End Build"

include $(KOKKOS_PATH)/Makefile.kokkos

# KOKKOS_INTERNAL_USE_CUDA is not exported to installed Makefile.kokkos
# use KOKKOS_DEVICE here
ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
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
	-I${SRC_DIR}/../common	\
	-I${SRC_DIR}/../fixture	\
	-I${SRC_DIR}


EXE_EXAMPLE_FENL = KokkosExample_Fenl
OBJ_EXAMPLE_FENL = BoxElemPart.o main.o fenl.o

TARGETS = $(EXE_EXAMPLE_FENL)

#TEST_TARGETS =

$(EXE_EXAMPLE_FENL) : $(OBJ_EXAMPLE_FENL) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_EXAMPLE_FENL) $(KOKKOS_LIBS) $(LIB) -o $(EXE_EXAMPLE_FENL)

build_all : $(TARGETS)


test : build_all

clean:
	rm -f *.o $(EXE_EXAMPLE_FENL) KokkosCore_config.*

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(EXAMPLE_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

