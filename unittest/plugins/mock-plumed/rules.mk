# Make rules of the PLUMED stand-in, included by the Makefile that "configure" writes.
# Where the files are put: src/lib/install/ of the build folder, as PLUMED does it (the
# Windows cross-compilation of LAMMPS takes them from there), and on "make install" the
# folders of an installed PLUMED.  The kernel (libplumedKernel) is the same code as the
# wrapper: nothing in it does any work.

INSTALL = src/lib/install
TARGETS = $(INSTALL)/libplumed.a $(INSTALL)/libplumedWrapper.a $(INSTALL)/plumed$(EXEEXT)
ifeq ($(SHARED),yes)
TARGETS += $(INSTALL)/libplumed.so $(INSTALL)/libplumedKernel.so
endif
INCLUDES = -I$(srcdir)/src/include

all: $(TARGETS)

$(INSTALL)/wrapper.o: $(srcdir)/src/lib/Plumed.cpp $(srcdir)/src/include/plumed/wrapper/Plumed.h
	mkdir -p $(INSTALL)
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) -c -o $@ $<

$(INSTALL)/wrapper-pic.o: $(srcdir)/src/lib/Plumed.cpp $(srcdir)/src/include/plumed/wrapper/Plumed.h
	mkdir -p $(INSTALL)
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) -fPIC -c -o $@ $<

$(INSTALL)/libplumed.a $(INSTALL)/libplumedWrapper.a: $(INSTALL)/wrapper.o
	rm -f $@
	$(AR) rc $@ $^
	$(RANLIB) $@

$(INSTALL)/libplumed.so $(INSTALL)/libplumedKernel.so: $(INSTALL)/wrapper-pic.o
	$(CXX) -shared -Wl,-soname,$(notdir $@) $(LDFLAGS) -o $@ $^

$(INSTALL)/plumed$(EXEEXT): $(srcdir)/src/lib/main.cpp
	mkdir -p $(INSTALL)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

install: all
	mkdir -p $(DESTDIR)$(libdir) $(DESTDIR)$(bindir) $(DESTDIR)$(includedir)/plumed/wrapper
	cp $(INSTALL)/lib*.* $(DESTDIR)$(libdir)/
	cp $(INSTALL)/plumed$(EXEEXT) $(DESTDIR)$(bindir)/
	cp $(srcdir)/src/include/plumed/wrapper/Plumed.h $(DESTDIR)$(includedir)/plumed/wrapper/

.PHONY: all install
