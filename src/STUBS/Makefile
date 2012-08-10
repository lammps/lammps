# Makefile for MPI stubs library

# Syntax:
#   make                 # build static lib as libmpi_stubs.a
#   make shlib           # build shared lib as libmpi_stubs.so
#   make clean           # remove *.o and lib files

# edit System-specific settings as needed for your platform

SHELL = /bin/sh
.IGNORE:

# Files

SRC =		mpi.c
INC =		mpi.h

# Definitions

EXE =		libmpi_stubs.a
SHLIB =		libmpi_stubs.so
OBJ = 		$(SRC:.c=.o)

# System-specific settings

CC =		g++
CCFLAGS =	-O
SHFLAGS =	-fPIC

ARCHIVE =	ar
ARCHFLAG =	rs
SHLIBFLAGS =    -shared          

# Targets

lib:	$(OBJ)
	$(ARCHIVE) $(ARCHFLAG) $(EXE) $(OBJ)

shlib:	$(OBJ)
	$(CC) $(CFLAGS) $(SHFLAGS) $(SHLIBFLAGS) -o $(SHLIB) $(OBJ)

clean:
	rm -f *.o libmpi_stubs.a libmpi_stubs.so

# Compilation rules

.c.o:
	$(CC) $(CCFLAGS) $(SHFLAGS) -c $<

# Individual dependencies

$(OBJ):	$(INC)
