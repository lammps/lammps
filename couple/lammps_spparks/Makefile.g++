# Makefile for grain growth via LAMMPS <-> SPPARKS coupling

SHELL = /bin/sh

# System-specific settings

LAMMPS =	/home/sjplimp/lammps
SPPARKS =	/home/sjplimp/spparks

CC =		g++
CCFLAGS =	-g -O -DMPICH_IGNORE_CXX_SEEK -I../library
DEPFLAGS =	-M
LINK =		g++
LINKFLAGS =	-g -O -L../library -L${LAMMPS}/src -L${SPPARKS}/src
USRLIB =	-lcouple -llmp_g++ -lspk_g++ 
SYSLIB =	-lfftw -lmpich -lpthread
ARCHIVE =	ar
ARFLAGS =	-rc
SIZE =		size

# Files

EXE = 	lmpspk
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# Targets

lmpspk:	$(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(USRLIB) $(SYSLIB) -o $(EXE)
	$(SIZE) $(EXE)

clean:
	rm $(EXE) *.o

# Compilation rules

%.o:%.cpp
	$(CC) $(CCFLAGS) -c $<

%.d:%.cpp
	$(CC) $(CCFLAGS) $(DEPFLAGS) $< > $@

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
include $(DEPENDS)
