KOKKOS_DEVICES=OpenMP
KOKKOS_CUDA_OPTIONS=enable_lambda
KOKKOS_ARCH = "SNB,Kepler35"


MAKEFILE_PATH := $(subst Makefile,,$(abspath $(lastword $(MAKEFILE_LIST))))

ifndef KOKKOS_PATH
  KOKKOS_PATH = $(MAKEFILE_PATH)../..
endif

SRC = $(wildcard $(MAKEFILE_PATH)*.cpp)
HEADERS = $(wildcard $(MAKEFILE_PATH)*.hpp)

vpath %.cpp $(sort $(dir $(SRC)))

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = make_buildlink.cuda
else
CXX = g++
EXE = make_buildlink.host
endif

CXXFLAGS ?= -O3 -g
override CXXFLAGS += -I$(MAKEFILE_PATH)

DEPFLAGS = -M
LINK = ${CXX}
LINKFLAGS =

OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

test: build
	./$(EXE)

$(EXE): $(SRC) $(KOKKOS_LINK_DEPENDS) $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(KOKKOS_CXXLDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(SRC) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean
	rm -f *.o *.cuda *.host

