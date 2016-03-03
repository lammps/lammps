# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Tell make that these are phony targets
.PHONY: all help clean install uninstall

include config.mk

# Build all of the executable files
all:
	$(MAKE) -C src
	$(MAKE) -C examples

# Build the help files (with Doxygen)
help:
	$(MAKE) -C src help

# Clean up the executable files
clean:
	$(MAKE) -C src clean
	$(MAKE) -C examples clean

# Install the executable, man page, and shared library
install:
	$(MAKE) -C src
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/bin
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/lib
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/man
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/man/man1
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/include
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS_EXEC) src/voro++ $(PREFIX)/bin
	$(INSTALL) $(IFLAGS) man/voro++.1 $(PREFIX)/man/man1
	$(INSTALL) $(IFLAGS) src/libvoro++.a $(PREFIX)/lib
	$(INSTALL) $(IFLAGS) src/voro++.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/c_loops.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/cell.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/common.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/config.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/container.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/container_prd.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/rad_option.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/pre_container.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/unitcell.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/v_base.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/v_compute.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/wall.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/worklist.hh $(PREFIX)/include/voro++

# Uninstall the executable, man page, and shared library
uninstall:
	rm -f $(PREFIX)/bin/voro++
	rm -f $(PREFIX)/man/man1/voro++.1
	rm -f $(PREFIX)/lib/libvoro++.a
	rm -f $(PREFIX)/include/voro++/voro++.hh
	rm -f $(PREFIX)/include/voro++/c_loops.hh
	rm -f $(PREFIX)/include/voro++/cell.hh
	rm -f $(PREFIX)/include/voro++/common.hh
	rm -f $(PREFIX)/include/voro++/config.hh
	rm -f $(PREFIX)/include/voro++/container.hh
	rm -f $(PREFIX)/include/voro++/container_prd.hh
	rm -f $(PREFIX)/include/voro++/pre_container.hh
	rm -f $(PREFIX)/include/voro++/rad_option.hh
	rm -f $(PREFIX)/include/voro++/unitcell.hh
	rm -f $(PREFIX)/include/voro++/v_base.hh
	rm -f $(PREFIX)/include/voro++/v_compute.hh
	rm -f $(PREFIX)/include/voro++/wall.hh
	rm -f $(PREFIX)/include/voro++/worklist.hh
	rmdir $(PREFIX)/include/voro++
