SHELL = /bin/sh

# Path to LAMMPS extraction directory
LAMMPS_ROOT = ../svn-dist
LAMMPS_SRC = $(LAMMPS_ROOT)/src

# Remove the line below if using mpicxx/mpic++ as your C++ compiler
MPI_STUBS = $(LAMMPS_SRC)/STUBS

FC = gfortran # replace with your Fortran compiler
CXX = g++     # replace with your C++ compiler

# Flags for Fortran compiler, C++ compiler, and C preprocessor, respectively
FFLAGS = -O2
CXXFLAGS = -O2
CPPFLAGS =

all : liblammps_fortran.a

liblammps_fortran.a : LAMMPS.o LAMMPS-wrapper.o
	$(AR) rs $@ $^

LAMMPS.o lammps.mod : LAMMPS.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $<

LAMMPS-wrapper.o : LAMMPS-wrapper.cpp LAMMPS-wrapper.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -I$(LAMMPS_SRC) -I$(MPI_STUBS)

clean :
	$(RM) *.o *.mod liblammps_fortran.a

dist :
	tar -czf Fortran-interface.tar.gz LAMMPS-wrapper.h LAMMPS-wrapper.cpp LAMMPS.F90 makefile README
