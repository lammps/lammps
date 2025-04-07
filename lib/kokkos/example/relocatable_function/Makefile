KOKKOS_PATH = ../..
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
SRC = $(wildcard ${KOKKOS_SRC_PATH}/example/relocatable_function/*.cpp)
vpath %.cpp $(sort $(dir $(SRC)))

KOKKOS_USE_DEPRECATED_MAKEFILES=1

default: build
	echo "Start Build"

EXE = Kokkos_example_relocatable_function

LINK = ${CXX}

DEPFLAGS = -M

OBJ = $(notdir $(SRC:.cpp=.o))

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

test: $(EXE)
	./$(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) -o $(EXE)

clean: kokkos-clean 
	rm -f *.o Kokkos_example_relocatable_function

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)
