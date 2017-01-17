KOKKOS_PATH ?= ../..

MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
SRC_DIR := $(dir $(MAKEFILE_PATH))

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=%.o)

#SRC = $(wildcard *.cpp)
#OBJ = $(SRC:%.cpp=%.o)

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
  CXX = $(KOKKOS_PATH)/bin/nvcc_wrapper
  EXE = $(addsuffix .cuda, $(shell basename $(SRC_DIR)))
else
  CXX = g++
  EXE = $(addsuffix .host, $(shell basename $(SRC_DIR)))
endif

CXXFLAGS = -O3 -I$(SRC_DIR)
LINK ?= $(CXX)
LDFLAGS ?=

include $(KOKKOS_PATH)/Makefile.kokkos

DEPFLAGS = -M

LIB =


build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: 
	rm -f *.a *.o *.cuda *.host

# Compilation rules

%.o:$(SRC_DIR)/%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

