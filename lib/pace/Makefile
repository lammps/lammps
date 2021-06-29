SHELL = /bin/sh

# ------ FILES ------

SRC_FILES = $(wildcard src/USER-PACE/*.cpp)
SRC = $(filter-out src/USER-PACE/pair_pace.cpp, $(SRC_FILES))

# ------ DEFINITIONS ------

LIB = libpace.a
OBJ =   $(SRC:.cpp=.o)


# ------ SETTINGS ------
CXXFLAGS = -O3 -fPIC -Isrc/USER-PACE

ARCHIVE =	ar
ARCHFLAG =	-rc
USRLIB =
SYSLIB =

# ------ MAKE PROCEDURE ------

lib: $(OBJ)
	$(ARCHIVE) $(ARFLAGS) $(LIB) $(OBJ)

# ------ COMPILE RULES ------

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ------ CLEAN ------
clean-all:
	-rm -f *~ $(OBJ) $(LIB)

clean-build:
	-rm -f *~ $(OBJ)
