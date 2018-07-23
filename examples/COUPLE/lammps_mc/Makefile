# Makefile for MC

SHELL = /bin/sh

SRC =	mc.cpp random_park.cpp
OBJ = 	$(SRC:.cpp=.o)

# change this line for your machine to path for CSlib src dir

CSLIB = /home/sjplimp/lammps/lib/message/cslib/src    

# compiler/linker settings

CC =		g++
CCFLAGS =	-g -O3 -I$(CSLIB)
LINK =		g++
LINKFLAGS =	-g -O -L$(CSLIB)

# targets

mc:	$(OBJ)
#	this line if built the CSlib within lib/message with ZMQ support
#       note this is using the serial (no-mpi) version of the CSlib
	$(LINK) $(LINKFLAGS) $(OBJ) -lcsnompi -lzmq -o mc
#	this line if built the CSlib without ZMQ support
#	$(LINK) $(LINKFLAGS) $(OBJ) -lcsnompi -o mc

clean:
	@rm -f *.o mc

# rules

%.o:%.cpp
	$(CC) $(CCFLAGS) -c $<
