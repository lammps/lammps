# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=loops polygons odd_even find_voro_cell

# Makefile rules
all: $(EXECUTABLES)

loops: loops.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o loops loops.cc -lvoro++

polygons: polygons.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o polygons polygons.cc -lvoro++

odd_even: odd_even.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o odd_even odd_even.cc -lvoro++

find_voro_cell: find_voro_cell.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o find_voro_cell find_voro_cell.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
