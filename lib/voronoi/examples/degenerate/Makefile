# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=degenerate degenerate2

# Makefile rules
all: degenerate degenerate2

degenerate: $(SOURCE) degenerate.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o degenerate degenerate.cc -lvoro++

degenerate2: $(SOURCE) degenerate2.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o degenerate2 degenerate2.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
