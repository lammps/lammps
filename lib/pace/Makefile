SHELL = /bin/sh

# ------ FILES ------

YAML_CPP_PATH = src/yaml-cpp
YAML_CPP_INC = $(YAML_CPP_PATH)/include

WIGNER_CPP_INC = src/wigner-cpp/include

CNPY_CPP_PATH = src/cnpy
CNPY_CPP_INC = src
CNPY_SRC_FILES =  $(CNPY_CPP_PATH)/cnpy.cpp

SRC_FILES = $(wildcard src/ML-PACE/ace/*.cpp) $(wildcard src/ML-PACE/ace-evaluator/*.cpp)
SRC = $(filter-out src/ML-PACE/pair_pace.cpp, $(SRC_FILES)) $(CNPY_SRC_FILES)

# ------ DEFINITIONS ------

LIB = libpace.a
OBJ =   $(SRC:.cpp=.o)


# ------ SETTINGS ------
CXXFLAGS = -O3 -fPIC -Isrc/ML-PACE -I$(YAML_CPP_INC) -I$(WIGNER_CPP_INC) -I$(CNPY_CPP_INC) -DEXTRA_C_PROJECTIONS -DCOMPUTE_B_GRAD

ARCHIVE =	ar
ARCHFLAG =	-rc
USRLIB =
SYSLIB =

# ------ MAKE PROCEDURE ------

lib: $(OBJ) lib-yaml-cpp
	$(ARCHIVE) $(ARFLAGS) $(LIB) $(OBJ)


lib-yaml-cpp:
	cd $(YAML_CPP_PATH) && $(MAKE) lib

# ------ COMPILE RULES ------

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ------ CLEAN ------
clean-all:
	-rm -f *~ $(OBJ) $(LIB)
	cd $(YAML_CPP_PATH) && $(MAKE) clean-all

clean-build:
	-rm -f *~ $(OBJ)
	cd $(YAML_CPP_PATH) && $(MAKE) clean-build

