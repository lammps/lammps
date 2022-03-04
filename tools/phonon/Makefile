.SUFFIXES : .o .cpp
# compiler and flags
CC     = g++ -Wno-unused-result
LINK   = $(CC) -static
CFLAGS = -O3 $(DEBUG) $(UFLAG)

OFLAGS = -O3 $(DEBUG)
INC    = $(LPKINC) $(TCINC) $(SPGINC) $(FFTINC)
LIB    = $(LPKLIB) $(TCLIB) $(SPGLIB) $(FFTLIB)

# cLapack library needed
LPKINC = -I/opt/clapack/3.2.1/include
LPKLIB = -L/opt/clapack/3.2.1/lib -lclapack -lblas -lf2c -lm

# Tricubic library needed
TCINC = -I/opt/tricubic/1.0/include
TCLIB = -L/opt/tricubic/1.0/lib -ltricubic

# spglib, used to get the irreducible q-points
# if SFLAG is not set, spglib won't be used.
SFLAG  = -DUseSPG
SPGINC = -I/opt/spglib/1.9.7/include/spglib
SPGLIB = -L/opt/spglib/1.9.7/lib -lsymspg

# FFTW 3， used to deduce the force constants in real space
# if FFLAG is not set, fftw won't be used.
FFLAG  = -DFFTW3
FFTINC = -I/opt/fftw/3.3.7/include
FFTLIB = -L/opt/fftw/3.3.7/lib -lfftw3

# Debug flags
# DEBUG = -g -DDEBUG
UFLAG = $(SFLAG) $(FFLAG)

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
	rm -f ${ROOT}.tar.gz; tar -czvf ${ROOT}.tar.gz *.cpp  *.h Makefile README

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
