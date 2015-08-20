KOKKOS_PATH = ../..

vpath %.cpp ${KOKKOS_PATH}/example/fixture ${KOKKOS_PATH}/example/fenl

EXAMPLE_HEADERS = $(wildcard $(KOKKOS_PATH)/example/common/*.hpp ${KOKKOS_PATH}/example/fixture/*.hpp ${KOKKOS_PATH}/example/fenl/*.hpp)

default: build_all
	echo "End Build"
        
include $(KOKKOS_PATH)/Makefile.kokkos

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
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
	-I${KOKKOS_PATH}/example/common	\
	-I${KOKKOS_PATH}/example/fixture	\
	-I${KOKKOS_PATH}/example/fenl


EXE_EXAMPLE_FENL = KokkosExample_Fenl
OBJ_EXAMPLE_FENL = BoxElemPart.o main.o fenl.o

TARGETS = $(EXE_EXAMPLE_FENL)

#TEST_TARGETS =

$(EXE_EXAMPLE_FENL) : $(OBJ_EXAMPLE_FENL) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_EXAMPLE_FENL) $(KOKKOS_LIBS) $(LIB) -o $(EXE_EXAMPLE_FENL)

build_all : $(TARGETS)


test : build_all


# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(EXAMPLE_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

