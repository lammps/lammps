KOKKOS_PATH ?= ../..

MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
SRC_DIR := $(dir $(MAKEFILE_PATH))

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=%.o)

#SRC = $(wildcard *.cpp)
#OBJ = $(SRC:%.cpp=%.o)

default: build
	echo "Start Build"

# use installed Makefile.kokkos
include $(KOKKOS_PATH)/Makefile.kokkos

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = $(NVCC_WRAPPER)
CXXFLAGS = -I$(SRC_DIR) -I$(CUDA_PATH) -O3
LINK = $(CXX)
LINKFLAGS = -L$(CUDA_PATH)/lib64 -lcusparse
EXE = $(addsuffix .cuda, $(shell basename $(SRC_DIR)))
#KOKKOS_DEVICES = "Cuda,OpenMP"
#KOKKOS_ARCH = "SNB,Kepler35"
else
CXX = g++
CXXFLAGS = -I$(SRC_DIR) -O3
LINK = $(CXX)
LINKFLAGS =  
EXE = $(addsuffix .host, $(shell basename $(SRC_DIR)))
#KOKKOS_DEVICES = "OpenMP"
#KOKKOS_ARCH = "SNB"
endif

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

