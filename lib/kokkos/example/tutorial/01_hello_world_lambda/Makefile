KOKKOS_PATH = ../../..
SRC = $(wildcard *.cpp)

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = nvcc_wrapper
CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS = 
EXE = $(SRC:.cpp=.cuda)
KOKKOS_DEVICES = "Cuda,OpenMP"
KOKKOS_ARCH = "SNB,Kepler35"
KOKKOS_CUDA_OPTIONS = "enable_lambda"
else
CXX = g++
CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS =  
EXE = $(SRC:.cpp=.host)
KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "SNB"
endif

DEPFLAGS = -M

OBJ = $(SRC:.cpp=.o)
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean 
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<
