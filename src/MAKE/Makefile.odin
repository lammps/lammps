# odin = 1400 cluster, g++, MPICH, no FFTs

SHELL = /bin/sh
#.IGNORE:

# System-specific settings

CC =		g++
CCFLAGS =	-O -I/opt/mpich-mx/include -DFFT_NONE -DGZIP
DEPFLAGS =	-M
LINK =		g++
LINKFLAGS =	-O -L/opt/mpich-mx/lib -L/opt/mx/lib
USRLIB =	-lmpich -lmyriexpress
SYSLIB =
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
