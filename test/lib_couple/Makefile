# Hey emacs, this is a -*- Makefile -*- !
SHELL=/bin/sh
CSRC=$(wildcard *.c)
XSRC=$(wildcard *.cpp)
EXE=$(CSRC:%.c=%.x) $(XSRC:%.cpp=%.x)
OUTPUTS=$(CSRC:%.c=log.%-$(MACH)$(TAG)) $(XSRC:%.cpp=log.%-$(MACH)$(TAG))
LIB=../../src/liblammps_$(MACH).a

default:
	@echo Run tests with "make test COMP=<MPI C/C++ compiler> MACH=<build make target> MPICMD=<MPI launch command> LMPFLAGS=<additional flags for LAMMPS>"

test: $(LIB) $(OUTPUTS)

log.%-$(MACH)$(TAG): %.x
	$(MPICMD) python  $< $(LMPFLAGS) -log none -screen $@ -echo screen 
	test -f log.lammps && mv log.lammps $@ || true

%.x	: %.c $(LIB)
	$(COMP) -I../../src $< $(LIB) -o $@

%.x	: %.cpp $(LIB)
	$(COMP) -I../../src $< $(LIB) -o $@

clean:
	-rm -f log.* lammps.py *.pyc liblammps*.so
	
