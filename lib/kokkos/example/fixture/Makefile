KOKKOS_PATH = ../..
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
vpath %.cpp ${KOKKOS_SRC_PATH}/example/fixture

EXAMPLE_HEADERS = $(wildcard $(KOKKOS_SRC_PATH)/example/common/*.hpp ${KOKKOS_SRC_PATH}/example/fixture/*.hpp )

default: build_all
	echo "End Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
  CXX = $(KOKKOS_PATH)/bin/nvcc_wrapper
else
  CXX = g++
endif

CXXFLAGS = -O3
LINK ?= $(CXX)
LDFLAGS ?=

include $(KOKKOS_PATH)/Makefile.kokkos        

KOKKOS_CXXFLAGS +=	\
	-I${KOKKOS_SRC_PATH}/example/common	\
	-I${KOKKOS_SRC_PATH}/example/fixture

EXE_EXAMPLE_FIXTURE = KokkosExample_Fixture
OBJ_EXAMPLE_FIXTURE = Main.o TestFixture.o BoxElemPart.o

TARGETS = $(EXE_EXAMPLE_FIXTURE)

#TEST_TARGETS =

$(EXE_EXAMPLE_FIXTURE) : $(OBJ_EXAMPLE_FIXTURE) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ_EXAMPLE_FIXTURE) $(KOKKOS_LIBS) $(LIB) -o $(EXE_EXAMPLE_FIXTURE)

build_all : $(TARGETS)

test : build_all

clean: kokkos-clean
	rm -f *.o $(TARGETS)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(EXAMPLE_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<
