SHELL=/bin/sh

CXX=g++ -std=c++11
#CXX=g++
CXXFLAGS=-O -g -Wall -I../../src $(shell pkg-config --cflags readline)
#CXXFLAGS=-O -g -Wall -I../../src
LDFLAGS= -L../../src -Wl,-rpath,../../src -llammps $(shell pkg-config --libs readline)
#LDFLAGS= -L../../src -llammps -lreadline

lammps-shell: lammps-shell.o
	$(CXX) -o $@ $^ $(LDFLAGS)

lammps-shell.o: lammps-shell.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean:
	@rm -f lammps-shell lammps-shell.o core *~ .lammps-history log.lammps log.cite
