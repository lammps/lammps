#---------------------------------------------------------------------------------------------------
#
# This is Makefile for builing the executable TMDGen
#
# Alexey N. Volkov, University of Alabama, avolkov1@ua.edu, 2020, Version 13.00
#
#---------------------------------------------------------------------------------------------------
EXEPATH = .

F90 = ifort
F90FLAGS = -O3 -ipo
LDFLAGS =

OBJS = TPMLib.o TPMGeom.o TMDGenData.o TMDGen3D.o TMDGen.o
EXE = $(EXEPATH)/TMDGen

# compile and load
default:
	@echo " "
	@echo "Compiling Code of Program TMDGen"
	@echo "FORTRAN 90"
	$(MAKE) $(EXE)

$(EXE):	$(OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $(EXE)  $(OBJS)

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(F90FLAGS) -c $*.f90

clean:
	rm -f *.o

