# Makefile for CSlib = client/server messaging library
# type "make help" for options

SHELL = /bin/sh

# ----------------------------------------
# should only need to change this section
# compiler/linker settings
# ----------------------------------------

CC =		g++
CCFLAGS =	-g -O3 -DZMQ_$(ZMQ) -DMPI_$(MPI)
SHFLAGS =	-fPIC
ARCHIVE =       ar
ARCHFLAGS =     -rc
SHLIBFLAGS =	-shared

# files

LIB =	libcsmpi.a
SHLIB =	libcsmpi.so
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# build with ZMQ support or not

zmq = 		yes
ZMQ =		$(shell echo $(zmq) | tr a-z A-Z)

ifeq ($(ZMQ),YES)
  ZMQLIB = -lzmq
else
  CCFLAGS += -I./STUBS_ZMQ
endif

# build with MPI support or not

mpi = 		yes
MPI =		$(shell echo $(mpi) | tr a-z A-Z)

ifeq ($(MPI),YES)
  CC = mpicxx
else
  CCFLAGS += -I./STUBS_MPI
  LIB = libcsnompi.a
  SHLIB = libcsnompi.so
endif

# targets

shlib:	shlib_parallel shlib_serial

lib:	lib_parallel lib_serial

all:	shlib lib

help:
	@echo 'make                   default = shlib'
	@echo 'make shlib             build 2 shared CSlibs: parallel & serial'
	@echo 'make lib               build 2 static CSlibs: parallel & serial'
	@echo 'make all               build 4 CSlibs: shlib and lib'
	@echo 'make shlib_parallel    build shared parallel CSlib'
	@echo 'make shlib_serial      build shared serial CSlib'
	@echo 'make lib_parallel      build static parallel CSlib'
	@echo 'make lib_serial        build static serial CSlib'
	@echo 'make ... zmq=no        build w/out ZMQ support'
	@echo 'make clean             remove all *.o files'
	@echo 'make clean-all         remove *.o and lib files'
	@echo 'make tar               create a tarball, 2 levels up'

shlib_parallel:
	$(MAKE) clean
	$(MAKE) shared zmq=$(zmq) mpi=yes

shlib_serial:
	$(MAKE) clean
	$(MAKE) shared zmq=$(zmq) mpi=no

lib_parallel:
	$(MAKE) clean
	$(MAKE) static zmq=$(zmq) mpi=yes

lib_serial:
	$(MAKE) clean
	$(MAKE) static zmq=$(zmq) mpi=no

static:	$(OBJ)
	$(ARCHIVE) $(ARCHFLAGS) $(LIB) $(OBJ)

shared:	$(OBJ)
	$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) -o $(SHLIB) $(OBJ) $(ZMQLIB)

clean:
	@rm -f *.o *.pyc

clean-all:
	@rm -f *.o *.pyc lib*.a lib*.so

tar:
	cd ../..; tar cvf cslib.tar cslib/README cslib/LICENSE \
		cslib/doc cslib/src cslib/test

# rules

%.o:%.cpp
	$(CC) $(CCFLAGS) $(SHFLAGS) -c $<
