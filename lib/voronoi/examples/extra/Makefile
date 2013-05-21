# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=box_cut cut_region superellipsoid irregular

# Makefile rules
all: $(EXECUTABLES)

box_cut: box_cut.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o box_cut box_cut.cc -lvoro++

cut_region: cut_region.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o cut_region cut_region.cc -lvoro++

superellipsoid: superellipsoid.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o superellipsoid superellipsoid.cc -lvoro++

irregular: irregular.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o irregular irregular.cc -lvoro++

finite_sys: finite_sys.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o finite_sys finite_sys.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
