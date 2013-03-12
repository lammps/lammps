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
LPKINC = -I/opt/clapack/3.2.1/include
LPKLIB = -L/opt/clapack/3.2.1/lib -lclapack -lblas -lf2c #-lm
#
# Tricubic library needed
TCINC = -I/opt/tricubic/1.0/include
TCLIB = -L/opt/tricubic/1.0/lib -ltricubic
#
# spglib 0.7.1, used to get the irreducible q-points
# if UFLAG is not set, spglib won't be used.
UFLAG  = -DUseSPG
SPGINC = -I/opt/spglib/0.7.1/include
SPGLIB = -L/opt/spglib/0.7.1/lib -lsymspg
# if spglib > 0.7.1 is used, please 
# 1) modify file phonon.cpp, instruction can be found by searching 0.7.1
# 2) uncomment the following two lines
#SPGINC = -I/opt/spglib/1.1.2/include
#SPGLIB = -L/opt/spglib/1.1.2/lib -lsymspg

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
	@echo "#define VERSION `svn info|grep '^Revision'|cut -d: -f2`" > version.h

#====================================================================
.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	$(CC) $(CFLAGS) $(INC) -c $<
