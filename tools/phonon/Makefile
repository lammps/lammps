.SUFFIXES : .o .cpp
# compiler and flags
CC     = g++ -Wall
LINK   = $(CC)
CFLAGS = -O3 $(DEBUG) $(UFLAG)
#
OFLAGS = -O3 $(DEBUG)
INC    = $(LPKINC) $(TCINC) $(SPGINC)
LIB    = $(LPKLIB) $(TCLIB) $(SPGLIB)
#
# cLapack library needed
LPKINC = 
LPKLIB =-llapack
#
#
# spglib 1.8.2, used to get the irreducible q-points
# if UFLAG is not set, spglib won't be used.

# UFLAG  = -DUseSPG
# SPGINC = -I/opt/libs/spglib/1.8.2/include
# SPGLIB = -L/opt/libs/spglib/1.8.2/lib -lsymspg

# if spglib other than version 1.8.2 is used, please 
# modify file phonon.cpp, instruction can be found by searching 1.8.2

# Debug flags
#DEBUG = -g -DDEBUG
#====================================================================
ROOT   = phana
# executable name
EXE    = $(ROOT)
#====================================================================
# source and rules
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

#====================================================================
all:  ${EXE}

${EXE}: $(OBJ)
	$(LINK) $(OFLAGS) $(OBJ) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod ${EXE}

tar:
	rm -f ${ROOT}.tar; tar -czvf ${ROOT}.tar.gz *.cpp  *.h Makefile README

ver:
	@echo "#define VERSION `git log|grep '^commit'|wc -l`" > version.h

#====================================================================
.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	$(CC) $(CFLAGS) $(INC) -c $<

#====================================================================
# dependencies
disp.o: disp.cpp phonon.h dynmat.h memory.h interpolate.h green.h timer.h \
 global.h
dynmat.o: dynmat.cpp dynmat.h memory.h interpolate.h version.h global.h
green.o: green.cpp green.h memory.h global.h
interpolate.o: interpolate.cpp interpolate.h memory.h global.h
main.o: main.cpp dynmat.h memory.h interpolate.h phonon.h
memory.o: memory.cpp memory.h
phonon.o: phonon.cpp phonon.h dynmat.h memory.h interpolate.h green.h \
 timer.h global.h
timer.o: timer.cpp timer.h
