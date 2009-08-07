# odin = 1400 cluster, g++, MPICH, no FFTs

SHELL = /bin/sh
#.IGNORE:

# System-specific settings
# LINKFORT & FORTLIB settings not needed if LAMMPS Fortran libs not used

include		Makefile.package

CC =		g++
CCFLAGS =	$(PKGINC) -O -I/opt/mpich-mx/include -DFFT_NONE -DLAMMPS_GZIP
DEPFLAGS =	-M
LINK =		g++
LINKFORT =	
LINKFLAGS =	$(PKGPATH) $(LINKFORT) -O -L/opt/mpich-mx/lib -L/opt/mx/lib
USRLIB =	$(PKGLIB) -lmpich -lmyriexpress
FORTLIB =	
SYSLIB =$(FORTLIB) 
ARCHIVE =	ar
ARFLAGS =	-rc
SIZE =		size

# Link target

$(EXE):	$(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(USRLIB) $(SYSLIB) -o $(EXE)
	$(SIZE) $(EXE)

# Library target

lib:	$(OBJ)
	$(ARCHIVE) $(ARFLAGS) $(EXE) $(OBJ)

# Compilation rules

%.o:%.cpp
	$(CC) $(CCFLAGS) -c $<

%.d:%.cpp
	$(CC) $(CCFLAGS) $(DEPFLAGS) $< > $@

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
include $(DEPENDS)
