# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=cell_statistics custom_output radical

# Makefile rules
all: $(EXECUTABLES)

cell_statistics: cell_statistics.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o cell_statistics cell_statistics.cc -lvoro++

custom_output: custom_output.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o custom_output custom_output.cc -lvoro++

radical: radical.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o radical radical.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
