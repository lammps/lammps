#Set your Kokkos path to something appropriate
KOKKOS_PATH = ${HOME}/git/kokkos-github-repo
#KOKKOS_DEVICES = "Cuda"
#KOKKOS_ARCH = "Pascal60"
#KOKKOS_CUDA_OPTIONS = enable_lambda
KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "Power8"

SRC = stream-kokkos.cc

default: build
	echo "Start Build"
	
CXXFLAGS = -O3
#CXX = ${HOME}/git/kokkos-github-repo/bin/nvcc_wrapper
CXX = g++

LINK = ${CXX}

LINKFLAGS =  
EXE = stream-kokkos

DEPFLAGS = -M

OBJ = $(SRC:.cc=.o)
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean 
	rm -f *.o $(EXE)

# Compilation rules

%.o:%.cc $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<
