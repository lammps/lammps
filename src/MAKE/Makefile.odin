# odin = 1400 cluster, g++, MPICH, no FFTs

SHELL = /bin/sh

# ---------------------------------------------------------------------
# compiler/linker settings
# specify flags and libraries needed for your compiler

CC =		g++
CCFLAGS =	-O
DEPFLAGS =	-M
LINK =		G++
LINKFLAGS =	-O
LIB =           
ARCHIVE =	ar
ARFLAGS =	-rc
SIZE =		size

# ---------------------------------------------------------------------
# LAMMPS-specific settings
# specify settings for LAMMPS features you will use

# LAMMPS ifdef options, see doc/Section_start.html

LMP_INC =	-DLAMMPS_GZIP

# MPI library, can be src/STUBS dummy lib
# INC = path for mpi.h, MPI compiler settings
# PATH = path for MPI library
# LIB = name of MPI library

MPI_INC =       -I/opt/mpich-mx/include
MPI_PATH =	-L/opt/mpich-mx/lib -L/opt/mx/lib
MPI_LIB =	-lmpich -lmyriexpress

# FFT library, can be -DFFT_NONE if not using PPPM from KSPACE package
# INC = -DFFT_FFTW, -DFFT_INTEL, -DFFT_NONE, etc, FFT compiler settings
# PATH = path for FFT library
# LIB = name of FFT library

FFT_INC =       -DFFT_NONE
FFT_PATH = 
FFT_LIB =	

# JPEG library, only needed if -DLAMMPS_JPEG listed with LMP_INC
# INC = path for jpeglib.h
# PATH = path for JPEG library
# LIB = name of JPEG library

JPG_INC =       
JPG_PATH = 	
JPG_LIB =	

# additional system libraries needed by LAMMPS package libraries
# these settings are IGNORED if the corresponding LAMMPS package
#   (e.g. gpu, meam) is NOT included in the LAMMPS build
# SYSLIB = names of libraries
# SYSPATH = paths of libraries

gpu_SYSLIB =       -lcudart
meam_SYSLIB =      -lifcore -lsvml -lompstub -limf
reax_SYSLIB =      -lifcore -lsvml -lompstub -limf
user-atc_SYSLIB =  -lblas -llapack

gpu_SYSPATH =      -L/usr/local/cuda/lib64
meam_SYSPATH =     -L/opt/intel/fce/10.0.023/lib
reax_SYSPATH =     -L/opt/intel/fce/10.0.023/lib
user-atc_SYSPATH = 	

# ---------------------------------------------------------------------
# build rules and dependencies
# no need to edit this section

include	Makefile.package

EXTRA_INC = $(LMP_INC) $(PKG_INC) $(MPI_INC) $(FFT_INC) $(JPG_INC)
EXTRA_PATH = $(PKG_PATH) $(MPI_PATH) $(FFT_PATH) $(JPG_PATH) $(PKG_SYSPATH)
EXTRA_LIB = $(PKG_LIB) $(MPI_LIB) $(FFT_LIB) $(JPG_LIB) $(PKG_SYSLIB)

# Link target

$(EXE):	$(OBJ)
	$(LINK) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(EXTRA_LIB) $(LIB) -o $(EXE)
	$(SIZE) $(EXE)

# Library target

lib:	$(OBJ)
	$(ARCHIVE) $(ARFLAGS) $(EXE) $(OBJ)

# Compilation rules

%.o:%.cpp
	$(CC) $(CCFLAGS) $(EXTRA_INC) -c $<

%.d:%.cpp
	$(CC) $(CCFLAGS) $(EXTRA_INC) $(DEPFLAGS) $< > $@

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
include $(DEPENDS)
