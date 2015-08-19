KOKKOS_PATH = ../..

vpath %.cpp ${KOKKOS_PATH}/example/fixture

EXAMPLE_HEADERS = $(wildcard $(KOKKOS_PATH)/example/common/*.hpp ${KOKKOS_PATH}/example/fixture/*.hpp )

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
	-I${KOKKOS_PATH}/example/fixture

OBJ_EXAMPLE_FIXTURE = Main.o TestFixture.o BoxElemPart.o
EXE_EXAMPLE_FIXTURE = KokkosExample_Fixture

TARGETS = $(EXE_EXAMPLE_FIXTURE)

#TEST_TARGETS =

$(EXE_EXAMPLE_FIXTURE) : $(OBJ_EXAMPLE_FIXTURE) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_EXAMPLE_FIXTURE) $(KOKKOS_LIBS) $(LIB) -o $(EXE_EXAMPLE_FIXTURE)


build_all : $(TARGETS)


test : build_all


# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(EXAMPLE_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

