# Hey emacs, this is a -*- Makefile -*- !
SHELL=/bin/sh
DIRS=$(wildcard [0-9]*)
EXE=../src/lmp_$(MACH)

default:
	@echo Run tests with "make test MACH=<build make target> MPICMD=<MPI launch command> LMPFLAGS=<additional flags for LAMMPS>"

test: $(EXE)
	for dir in $(DIRS); do cd $$dir ; $(MAKE) $(MFLAGS) MACH=$(MACH) MPICMD="$(MPICMD)" LMPFLAGS="$(LMPFLAGS)" test || exit 1; cd ..; done

clean:
	for dir in $(DIRS); do cd $$dir ; $(MAKE) $(MFLAGS) clean; cd ..; done
	
