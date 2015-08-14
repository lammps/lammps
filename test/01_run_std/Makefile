# Hey emacs, this is a -*- Makefile -*- !
SHELL=/bin/sh
INPUTS=$(wildcard in.*)
OUTPUTS=$(INPUTS:in.%=log.%-$(MACH)$(TAG))
EXE=../../src/lmp_$(MACH)

default:
	@echo Run tests with "make test MACH=<build make target> MPICMD=<MPI launch command> LMPFLAGS=<additional flags for LAMMPS>"

test: $(EXE) $(OUTPUTS)

log.%-$(MACH)$(TAG): in.% $(EXE)
	$(MPICMD) $(EXE) $(LMPFLAGS) -log none -screen $@ -echo screen -in $<

clean:
	-rm -f log.*
	
