.SUFFIXES : .o .cpp
# compiler and flags
CC     = g++ -Wno-unused-result
LINK   = $(CC) -static
CFLAGS = -O3 $(DEBUG) $(UFLAG)
#
OFLAGS = -O3 $(DEBUG)
INC    = $(LPKINC) $(TCINC) $(SPGINC)
LIB    = $(LPKLIB) $(TCLIB) $(SPGLIB)
#
# cLapack library needed
LPKINC = -I/opt/libs/clapack/3.2.1/include
LPKLIB = -L/opt/libs/clapack/3.2.1/lib -lclapack -lblas -lf2c #-lm
#
# Tricubic library needed
TCINC = -I/opt/libs/tricubic/1.0/include
TCLIB = -L/opt/libs/tricubic/1.0/lib -ltricubic
#
# spglib 1.8.2, used to get the irreducible q-points
# if UFLAG is not set, spglib won't be used.
UFLAG  = -DUseSPG
SPGINC = -I/opt/libs/spglib/1.8.2/include
SPGLIB = -L/opt/libs/spglib/1.8.2/lib -lsymspg
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
all:  ver ${EXE}

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
