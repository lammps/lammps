KOKKOS_PATH = ../../../..
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
SRC = $(wildcard ${KOKKOS_SRC_PATH}/example/tutorial/Hierarchical_Parallelism/02_nested_parallel_for/*.cpp)
vpath %.cpp $(sort $(dir $(SRC)))

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS = 
EXE = 02_nested_parallel_for.cuda
KOKKOS_DEVICES = "Cuda,OpenMP"
KOKKOS_ARCH = "SNB,Kepler35"
else
CXX = g++
CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS =  
EXE = 02_nested_parallel_for.host
KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "SNB"
endif
ifneq (,$(findstring ROCm,$(KOKKOS_DEVICES)))
CXX = /opt/rocm/hcc/bin/clang++
CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS =
EXE = 02_nested_parallel_for.rocm
KOKKOS_DEVICES = "ROCm"
KOKKOS_ARCH = "Fiji"
endif


DEPFLAGS = -M

OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

test: $(EXE)
	./$(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean 
	rm -f *.o *.cuda *.host *.rocm

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)
