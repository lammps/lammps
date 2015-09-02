# Hey emacs, this is a -*- Makefile -*- !
SHELL=/bin/sh
INPUTS=$(wildcard in.*.py)
OUTPUTS=$(INPUTS:in.%.py=log.%-$(MACH)$(TAG))
LIB=liblammps.so lammps.py

default:
	@echo Run tests with "make test MACH=<build make target> MPICMD=<MPI launch command> LMPFLAGS=<additional flags for LAMMPS>"

test: $(LIB) $(OUTPUTS)

log.%-$(MACH)$(TAG): in.%.py $(LIB)
	$(MPICMD) python  $< $(LMPFLAGS) -log none -screen $@ -echo screen 
	test -f log.lammps && mv log.lammps $@ || true

clean:
	-rm -f log.* lammps.py *.pyc liblammps*.so
	
