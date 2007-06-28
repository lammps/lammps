CC 	= gcc
LN 	= gcc
OPTIM	= -O2 -finline-functions -g
XLIBDIR	= -L/usr/X11R6/lib
XINCLUDE = 
DEFINES	= -DMISSINGDEFS -D_POSIX_SOURCE -DUSEPRIVATE -DINCL_FLOAT
WARNS	= -Wimplicit -Wunused -Wmissing-prototypes
CCOPTS 	= -c $(OPTIM) $(DEFINES) $(WARNS) $(XINCLUDE)
LNOPTS	= $(OPTIM) 

LIBS	= $(XLIBDIR) -lX11 -lXaw -lm

# try this version on a Mac
#LIBS	= $(XLIBDIR) -lX11 -lXaw -lXt -lXext -lXmu -lXpm

OBJS	= xmovie.o \
	control.o \
	scene.o \
	read.o \
	hpsort.o

EXES = xmovie dummy tobinary toascii

.c.o:
	$(CC) $(CCOPTS) $*.c

xmovie: $(OBJS)
	$(CC) $(CCOPTS) version.c
	$(LN) -o xmovie $(LNOPTS) version.o $(OBJS) $(LIBS)

xmovie.o:   xmovie.h resource.h xmovie.c
control.o:  xmovie.h control.c
scene.o:    xmovie.h scene.c
read.o:     xmovie.h read.c
hpsort.o:   hpsort.c

clean:
	rm -f *.o *~ $(EXES)

tar:
	tar cvf xmovie.tar Makefile *.c *.h
	compress xmovie.tar
	uuencode xmovie.tar.Z xmovie.tar.Z > xmovie.tar.Z.uu

clean_obj:
	rm -f *.o

all: $(EXES)

dummy: dummy.o
	$(LN) $(LNOPTS) -o dummy dummy.o -lm

dummy.o: dummy.c

toascii: toascii.o
	$(LN) $(LNOPTS) -o toascii toascii.o

toascii.o: toascii.c

tobinary: tobinary.o
	$(LN) $(LNOPTS) -o tobinary tobinary.o

tobinary.o: tobinary.c
