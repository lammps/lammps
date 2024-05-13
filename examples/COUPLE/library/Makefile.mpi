# Makefile for coupling library

SHELL = /bin/sh

# System-specific settings

CC =		mpicxx
CCFLAGS =	-g -O -DMPICH_IGNORE_CXX_SEEK  -DOMPI_SKIP_MPICXX=1 
DEPFLAGS =	-M
LINK =		mpicxx
LINKFLAGS =	-g -O
ARCHIVE =	ar
ARFLAGS =	-rc
SIZE =		size

# Files

LIB = 	libcouple.a
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# Targets

lib:	$(OBJ)
	$(ARCHIVE) $(ARFLAGS) $(LIB) $(OBJ)

clean:
	rm $(LIB) *.o *.d

# Compilation rules

%.o:%.cpp
	$(CC) $(CCFLAGS) -c $<

%.d:%.cpp
	$(CC) $(CCFLAGS) $(DEPFLAGS) $< > $@

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
include $(DEPENDS)
