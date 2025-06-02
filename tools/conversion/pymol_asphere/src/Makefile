#***************************************************************************
#                               Makefile
#                          -------------------
#
#  _________________________________________________________________________
#    Build for the Graphics Utilities
#  _________________________________________________________________________
#
#    begin                : Thu June 9 2005
#    copyright            : (C) 2003 by W. Michael Brown
#    email                : wmbrown@sandia.gov
# ***************************************************************************/

#Compiler type
#COMPILER = intel
COMPILER = gnu
#COMPILER = mpi
#COMPILER = mingw

#Locations of outside objects relative to a source directory
HOBJ_DIR = ../obj

BIN_DIR = ../bin

ALL_DIR = .
ALL_LIB = $(HOBJ_DIR)/liball.a

GRPHICS_DIR = .
GRPHICS_LIB = $(HOBJ_DIR)/libgraphics.a

GRID_DIR = .
GRID_LIB = $(HOBJ_DIR)/libgrid.a

MOL_DIR = .
MOL_LIB = $(HOBJ_DIR)/libmol.a

MATH_DIR = .
MATH_LIB = $(HOBJ_DIR)/libmath.a

EVERY_LIB = $(MOL_LIB) $(GRID_LIB) $(GRPHICS_LIB) $(MATH_LIB) $(ALL_LIB)

MOLSIM_DIR = .
MOLSIM_LIB = $(HOBJ_DIR)/molsim.o $(HOBJ_DIR)/dynmif.o $(HOBJ_DIR)/dynmifq.o

# Include directories
INC         = -I$(ALL_DIR) -I$(MOLSIM_DIR) -I$(MOL_DIR) -I$(MATH_DIR) -I$(GRID_DIR) -I$(GRPHICS_DIR)

ifeq ($(COMPILER),intel)
  CPP         = icpc      # C++ Compiler
  CC          = icc       # C compiler
  AR          = xiar      #ar
  DBUG        = -g -DDEBUG -DNANCHECK #-Wall #-ansi
  OPT         = -O2 -xP -ipo -no-prec-div -static
endif

ifeq ($(COMPILER),gnu)
  CPP         = g++ # C++ Compiler
  CC          = gcc # C compiler
  AR          = ar
  DBUG        = -O2 -g -Wall #-DDEBUG -DNANCHECK -Wall -pedantic #-ansi
  OPT         = #-O3 -g -Wall
endif

ifeq ($(COMPILER),mpi)
  CPP         = mpic++ -DMUSE_MPI      # C++ Compiler
  CC          = mpicc -DMUSE_MPI    # C compiler
  AR          = ar
  DBUG        = -g -DDEBUG -DNANCHECK #-Wall #-pedantic #-ansi
  OPT         = #-xN -O3 #-ipo -no-prec-div -static #-O3
endif

ifeq ($(COMPILER),mingw)
  CPP         = /cygdrive/c/MINGW/bin/g++ # C++ Compiler
  CC          = /cygdrive/c/MINGW/bin/gcc # C compiler
  AR          = /cygdrive/c/MINGW/bin/ar
  DBUG        = #-g -DDEBUG -DNANCHECK -Wall -pedantic #-ansi
  OPT         = -O3 -static
endif

# Large file support?
LFSC        = #-D_LARGEFILE_SOURCE `getconf LFS_CFLAGS`
LFSL        = #`getconf LFS_LDFLAGS` `getconf LFS_LIBS`

# GNU Scientific Library?
GSLC        = #-DUSEGSL -I/usr/local/include/
GSLL        = #-lgsl -lgslcblas

# GA Library?
LIBGAC      = #-DLIBGA -I/usr/local/include/
LIBGAL      = #-lga -L/usr/local/lib/ -Wl,--allow-multiple-definition

# Movie frame support?
MOVIE       = #-DMOVIE

# VTK ?
VTKH        = #-DUSEVTK -I/usr/local/include/vtk-5.0
VTKL        = #-lvtkWidgets -lvtkHybrid -lvtkVolumeRendering -lvtkRendering -lvtkIO -lvtkGenericFiltering -lvtkGraphics -lvtkImaging -l vtkFiltering -lvtkCommon -L/usr/X11R6/lib/ -lGL -lXt -lSM -lICE -lX11 -lXext -lpthread -ldl

CFLAGS      = $(OPT) $(MOVIE) $(DBUG) $(INC) $(GSLC) $(LIBGAC) $(VTKH) -c
LFLAGS      = $(OPT)
LLIBS       = $(GSLL) $(LIBGAL) $(LFSL) $(VTKL)

# Distribution Directories
DIST_BIN  = /home/wmbrown/distbin/
DIST_MAN    = /home/wmbrown/cpp/manpages/man1/
DIST_DOC  = /home/wmbrown/cpp/doc/

OBJ_DIR = $(HOBJ_DIR)

# Objects for this project
THIS_OBJ =  $(OBJ_DIR)/asphere_vis.o $(GRPHICS_LIB) $(MATH_LIB) $(ALL_LIB)
EXECS = $(BIN_DIR)/asphere_vis

all: paths $(EXECS)

paths: $(HOBJ_DIR) $(BIN_DIR)

$(HOBJ_DIR) $(BIN_DIR):
	mkdir $@

libraries:
	cd $(ALL_DIR); make; cd $(MATH_DIR); make; cd $(GRID_DIR); make; \
        cd $(GRPHICS_DIR); make;

ALL_OBJS = $(OBJ_DIR)/error.o $(OBJ_DIR)/commandline.o \
	   $(OBJ_DIR)/misc.o

$(OBJ_DIR)/error.o: error.h error.cpp
	$(CPP) $(CFLAGS) -o $@ error.cpp

$(OBJ_DIR)/commandline.o: commandline.h commandline.cpp
	$(CPP) $(CFLAGS) -o $@ commandline.cpp

$(OBJ_DIR)/misc.o: misc.h misc.cpp
	$(CPP) $(CFLAGS) -o $@ misc.cpp

$(ALL_LIB): $(ALL_OBJS)
	$(AR) -crusv $(ALL_LIB) $(ALL_OBJS)

GRPHICS_O   = $(OBJ_DIR)/colors.o $(OBJ_DIR)/glsurface.o

$(OBJ_DIR)/colors.o: colors.h colors.cpp
	$(CPP) $(CFLAGS) -o $@ colors.cpp

$(OBJ_DIR)/glsurface.o: glsurface.h glsurface.cpp
	$(CPP) $(CFLAGS) -o $@ glsurface.cpp

$(GRPHICS_LIB): $(GRPHICS_O)
	$(AR) -crusv $(GRPHICS_LIB) $(GRPHICS_O)

MATH_OBJS      = $(OBJ_DIR)/cartesian.o $(OBJ_DIR)/miscm.o \
                 $(OBJ_DIR)/spherical.o 

$(OBJ_DIR)/miscm.o: miscm.h miscm.cpp
	$(CPP) $(CFLAGS) -o $@ miscm.cpp

$(OBJ_DIR)/cartesian.o: cartesian.h cartesian.cpp
	$(CPP) $(CFLAGS) -o $@ cartesian.cpp

$(OBJ_DIR)/spherical.o: spherical.h spherical.cpp
	$(CPP) $(CFLAGS) -o $@ spherical.cpp

$(MATH_LIB): $(MATH_OBJS)
	$(AR) -crusv $(MATH_LIB) $(MATH_OBJS)

$(OBJ_DIR)/asphere_vis.o: asphere_vis.cpp 
	$(CPP) $(CFLAGS) -o $@ asphere_vis.cpp

$(BIN_DIR)/asphere_vis: $(THIS_OBJ)
	$(CPP) $(LFLAGS) -o $@ $(THIS_OBJ) $(LLIBS)
#
#  Documentation
#
manpages: all
	/bin/tcsh make_manpages.sh

#
# Create a .tar distribution file
#
dist: all manpages
	/bin/tcsh makedistribution.sh


#
#  INSTALL to Mike's Directories
install: all manpages
	/bin/cp $(EXECS) $(DIST_BIN); \
	/bin/cp ./manpages/*.1 $(DIST_MAN); \
	/bin/cp ./manpages/*.pdf $(DIST_DOC)


#
#  Remove objects, cores, etc.
#

clean:
	rm -rf $(EXECS) $(THIS_OBJ) $(ALL_OBJ) $(MATH_OBJ) $(GRPHICS_O)
	cd $(OBJ_DIR); rm -f *.o

veryclean: clean
	rm -rf *~ ./api ./manpages

cleanproject: clean
	cd $(ALL_DIR); make clean; cd $(MATH_DIR); make clean; \
  cd $(GRPHICS_DIR); make clean;
