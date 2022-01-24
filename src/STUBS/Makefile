# Makefile for MPI stubs library

# Syntax:
#   make                 # build lib as libmpi_stubs.a
#   make clean           # remove *.o and lib files

# edit System-specific settings as needed for your platform

SHELL = /bin/sh
.IGNORE:

# Files

SRC =		mpi.cpp
INC =		mpi.h

# Definitions

EXE =		libmpi_stubs.a
OBJ = 		$(SRC:.cpp=.o)

# System-specific settings

CC =		g++
CCFLAGS =	-O -fPIC -I. # add -I to insure mpi.h from this dir is included
ARCHIVE =	ar
ARCHFLAG =	rs

# Targets

$(EXE):	$(OBJ)
	$(ARCHIVE) $(ARCHFLAG) $(EXE) $(OBJ)

clean:
	rm -f *.o libmpi_stubs.a

# Compilation rules

.cpp.o:
	$(CC) $(CCFLAGS) -c $<

# Individual dependencies

$(OBJ):	$(INC) ../version.h
