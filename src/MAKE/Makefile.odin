# odin = 1400 cluster, g++, MPICH, no FFTs

SHELL = /bin/sh
#.IGNORE:

# System-specific settings

# LINKFORT/FORTLIB settings can be removed if not using meam or reax packages
# LINKBLAS/BLASLIB settings can be removed if not using user-atc package
# LINKGPU/GPULIB settings can be removed if not using gpu package

include		Makefile.package

CC =		g++
CCFLAGS =	$(PKGINC) -O -I/opt/mpich-mx/include -DFFT_NONE -DLAMMPS_GZIP
DEPFLAGS =	-M
LINK =		g++
LINKFORT =	
LINKBLAS =	
LINKGPU =	
LINKFLAGS =	$(PKGPATH) $(LINKFORT) $(LINKBLAS) $(LINKGPU)
USRLIB =	$(PKGLIB) -lmpich -lmyriexpress
FORTLIB =	
BLASLIB =	
GPULIB =	
SYSLIB =$(FORTLIB) $(BLASLIB) $(GPULIB)
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
