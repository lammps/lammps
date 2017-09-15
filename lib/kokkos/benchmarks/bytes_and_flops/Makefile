KOKKOS_PATH = ${HOME}/kokkos
SRC = $(wildcard *.cpp)
KOKKOS_DEVICES=Cuda
KOKKOS_CUDA_OPTIONS=enable_lambda

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = bytes_and_flops.cuda
KOKKOS_DEVICES = "Cuda,OpenMP"
KOKKOS_ARCH = "SNB,Kepler35"
else
CXX = g++
EXE = bytes_and_flops.host
KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "SNB"
endif

CXXFLAGS = -O3 -g

DEPFLAGS = -M
LINK = ${CXX}
LINKFLAGS =

OBJ = $(SRC:.cpp=.o)
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) bench.hpp bench_unroll_stride.hpp bench_stride.hpp
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<
