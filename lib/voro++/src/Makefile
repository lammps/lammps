# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../config.mk

# List of the common source files
objs=cell.o common.o container.o unitcell.o v_compute.o c_loops.o \
     v_base.o wall.o pre_container.o container_prd.o
src=$(patsubst %.o,%.cc,$(objs))

# Makefile rules
all: libvoro++.a voro++

depend:
	$(CXX) -MM $(src) >Makefile.dep

include Makefile.dep

libvoro++.a: $(objs)
	rm -f libvoro++.a
	ar rs libvoro++.a $^

voro++: libvoro++.a cmd_line.cc
	$(CXX) $(CFLAGS) -L. -o voro++ cmd_line.cc -lvoro++

%.o: %.cc
	$(CXX) $(CFLAGS) -c $<

help: Doxyfile $(SOURCE)
	doxygen Doxyfile

clean:
	rm -f $(objs) voro++ libvoro++.a

.PHONY: all help execs depend
