# Hey emacs, this is a -*- Makefile -*- !
SHELL=/bin/sh
CSRC=$(wildcard *.c)
XSRC=$(wildcard *.cpp)
EXE=$(CSRC:%.c=%.x) $(XSRC:%.cpp=%.x)
OUTPUTS=$(CSRC:%.c=log.%-$(MACH)$(TAG)) $(XSRC:%.cpp=log.%-$(MACH)$(TAG))
LIB=../../src/liblammps_$(MACH).a

include ../../src/Makefile.package.settings
include ../../src/Makefile.package

default:
	@echo Run tests with "make test COMP=<MPI C/C++ compiler> MACH=<build make target> MPICMD=<MPI launch command> LMPFLAGS=<additional flags for LAMMPS>"

test: $(LIB) $(EXE) $(OUTPUTS)

log.%-$(MACH)$(TAG): %.x
	$(MPICMD) env PATH=".:$(PATH)" $< 1 $(<:%.x=in.%) #$(LMPFLAGS) -log none -screen $@ -echo screen 
	test -f log.lammps && mv log.lammps $@ || true

%.x	: %.c $(LIB)
	$(COMP) $(PKG_INC) $(PKG_SYSINC) $(LMP_INC) -I../../src $< $(LIB) $(PKG_PATH) $(PKG_SYSPATH) $(PKG_LIB) $(PKG_SYSLIB) $(JPG_LIB) -o $@

%.x	: %.cpp $(LIB)
	$(COMP) $(PKG_INC) $(PKG_SYSINC) $(LMP_INC) -I../../src $< $(LIB) $(PKG_PATH) $(PKG_SYSPATH) $(PKG_LIB) $(PKG_SYSLIB) $(JPG_LIB) -o $@

clean:
	-rm -f log.* lammps.py *.pyc liblammps*.so *.x
	
