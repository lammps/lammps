# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=single_cell platonic random_points import

# Makefile rules
all: $(EXECUTABLES)

single_cell: single_cell.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o single_cell single_cell.cc -lvoro++

platonic: platonic.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o platonic platonic.cc -lvoro++

random_points: random_points.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o random_points random_points.cc -lvoro++

import: import.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o import import.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
