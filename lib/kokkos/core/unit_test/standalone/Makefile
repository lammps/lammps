KOKKOS_DEVICES=OpenMP
KOKKOS_CUDA_OPTIONS=enable_lambda
KOKKOS_ARCH = "SNB,Kepler35"

MAKEFILE_PATH := $(subst Makefile,,$(abspath $(lastword $(MAKEFILE_LIST))))

ifndef KOKKOS_PATH
  KOKKOS_PATH = $(MAKEFILE_PATH)../../../
endif

SRC = $(wildcard $(MAKEFILE_PATH)*.cpp)
SRC += $(MAKEFILE_PATH)/../TestStackTrace.cpp
HEADERS = $(wildcard $(MAKEFILE_PATH)*.hpp)
HEADERS = $(wildcard $(MAKEFILE_PATH)/../*.hpp)

vpath %.cpp $(sort $(dir $(SRC)))

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = test.cuda
else
CXX = g++
EXE = test.host
endif

CXXFLAGS ?= -O0 -g -rdynamic
override CXXFLAGS += -I$(MAKEFILE_PATH) -I$(KOKKOS_PATH)/core/unit_test -I$(KOKKOS_PATH)/tpls/gtest -DTESTFILE=$(TESTFILE)
#SRC += $(KOKKOS_PATH)/tpls/gtest/gtest/gtest-all.cc

DEPFLAGS = -M
LINK = ${CXX}
LINKFLAGS = -rdynamic

OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS) gtest-all.o
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) gtest-all.o  -o $(EXE)

clean: kokkos-clean
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)

gtest-all.o:$(KOKKOS_PATH)/tpls/gtest/gtest/gtest-all.cc
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $(KOKKOS_PATH)/tpls/gtest/gtest/gtest-all.cc
