# Hey emacs, this is a -*- Makefile -*- !
SHELL=/bin/sh
DIRS=$(wildcard [0-9]*)
EXE=../src/lmp_$(MACH)
LIB=../src/liblammps_$(MACH).so

default:
	@echo Run tests with "make test MACH=<build make target> MPICMD=<MPI launch command> LMPFLAGS=<additional flags for LAMMPS>"

test: $(EXE)
	for dir in $(DIRS); do cd $$dir ; $(MAKE) $(MFLAGS) MACH=$(MACH) MPICMD="$(MPICMD)" LMPFLAGS="$(LMPFLAGS)" test || exit 1; cd ..; done

python: lib_python/lammps.py lib_python/liblammps.so  $(LIB)
	for dir in lib_python; do cd $$dir ; $(MAKE) $(MFLAGS) MACH=$(MACH) MPICMD="$(MPICMD)" LMPFLAGS="$(LMPFLAGS)" test || exit 1; cd ..; done

lib_python/lammps.py lib_python/liblammps.so : ../python/install.py ../python/lammps.py $(LIB)
	python ../python/install.py lib_python

clean:
	for dir in $(DIRS) lib_python; do cd $$dir ; $(MAKE) $(MFLAGS) clean; cd ..; done
	
