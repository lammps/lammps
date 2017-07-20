SHELL = /bin/sh

# Path to LAMMPS extraction directory
LAMMPS_ROOT = ../../..
LAMMPS_SRC = $(LAMMPS_ROOT)/src

# Uncomment the line below if using the MPI stubs library
MPI_STUBS = #-I$(LAMMPS_SRC)/STUBS

FC = mpif90 # replace with your Fortran compiler
CXX = mpicc # replace with your C++ compiler

# Flags for Fortran compiler, C++ compiler, and C preprocessor, respectively
FFLAGS = -O2 -fPIC 
CXXFLAGS = -O2 -fPIC
CPPFLAGS = -DOMPI_SKIP_MPICXX=1 -DMPICH_SKIP_MPICXX

all : liblammps_fortran.a liblammps_fortran.so simpleF.x

liblammps_fortran.so : LAMMPS.o LAMMPS-wrapper.o LAMMPS-wrapper2.o
	$(FC) $(FFLAGS) -shared -o $@ $^

simpleF.x: simple.o LAMMPS.o LAMMPS-wrapper.o LAMMPS-wrapper2.o
	$(FC) $(FFLAGS) simple.o -o simpleF.x liblammps_fortran.a $(LAMMPS_SRC)/liblammps_mvapich.a -lstdc++ /usr/lib64/libfftw3.a 

liblammps_fortran.a : LAMMPS.o LAMMPS-wrapper.o LAMMPS-wrapper2.o
	$(AR) rs $@ $^

LAMMPS.o lammps.mod : LAMMPS.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $<

simple.o : simple.f90
	$(FC) $(FFLAGS) -c $<

LAMMPS-wrapper.o : LAMMPS-wrapper.cpp LAMMPS-wrapper.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -I$(LAMMPS_SRC) $(MPI_STUBS)

LAMMPS-wrapper2.o : LAMMPS-wrapper2.cpp LAMMPS-wrapper2.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -I$(LAMMPS_SRC) $(MPI_STUBS)

clean :
	$(RM) *.o *.mod liblammps_fortran.a liblammps_fortran.so

dist :
	tar -czvf fortran-interface-callback.tar.gz LAMMPS-wrapper.h LAMMPS-wrapper.cpp LAMMPS-wrapper2.h LAMMPS-wrapper2.cpp LAMMPS.F90 makefile README simple.f90
