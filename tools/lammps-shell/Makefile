SHELL=/bin/sh

CXX=g++
CXXFLAGS=-O -g -Wall -I../../src -DLAMMPS_LIB_NO_MPI
LDFLAGS= -L../../src -llammps -lreadline

lammps-shell: lammps-shell.o
	$(CXX) -o $@ $^ $(LDFLAGS)

lammps-shell.o: lammps-shell.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean:
	@rm -f lammps-shell lammps-shell.o core *~ .lammps-history
