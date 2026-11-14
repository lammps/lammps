# Make rules of the MBX stand-in, included by the Makefile that "configure" writes.
# "make install" honors DESTDIR and the folder variables, as the Windows cross-compilation
# of LAMMPS sets them.

INCLUDES = -I$(srcdir)/src
ifeq ($(MPI),yes)
CPPFLAGS += -DHAVE_MPI=1
endif

all: libmbx.a

system.o: $(srcdir)/src/bblock/system.cpp $(srcdir)/src/bblock/system.h
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) -c -o $@ $<

libmbx.a: system.o
	rm -f $@
	$(AR) rc $@ $^
	$(RANLIB) $@

install: all
	mkdir -p $(DESTDIR)$(libdir) $(DESTDIR)$(includedir)/bblock
	cp libmbx.a $(DESTDIR)$(libdir)/
	cp $(srcdir)/src/bblock/system.h $(DESTDIR)$(includedir)/bblock/

.PHONY: all install
