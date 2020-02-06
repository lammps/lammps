#---------------------------------------------------------------------------------------------------
#
# This is Makefile for builing the executable TMDPotGen
#
# Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, 2020, Version 13.00
#
#---------------------------------------------------------------------------------------------------
EXEPATH = .

F90 = ifort
F90FLAGS = -Ofast -mcmodel=medium
#F90 = pgf90
#F90FLAGS = -fast -mcmodel=medium
LDFLAGS =

OBJS = TPMLib.o LinFun2.o Spline1.o Spline2.o TPMGeom.o TubePotBase.o TubePotTrue.o TubePotMono.o TMDPotGen.o
EXE = $(EXEPATH)/TMDPotGen

# compile and load
default:
	@echo " "
	@echo "Compiling Code of Program TMDPotGen"
	@echo "FORTRAN 90"
	$(MAKE) $(EXE)

$(EXE):	$(OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $(EXE)  $(OBJS)

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(F90FLAGS) -c $*.f90

clean:
	rm -f *.o

