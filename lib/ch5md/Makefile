EXTRAMAKE=Makefile.lammps.empty

CC=h5cc
CFLAGS=-D_DEFAULT_SOURCE -O2
HDF5_PATH=/usr
INC=-I include
AR=ar
ARFLAGS=rc
LIB=libch5md.a

all: lib Makefile.lammps

build:
	mkdir -p build
build/ch5md.o: src/ch5md.c | build
	$(CC) $(INC) $(CFLAGS) -c $< -o $@

Makefile.lammps:
	cp Makefile.lammps.empty $@

.PHONY: all lib clean

$(LIB): build/ch5md.o
	$(AR) $(ARFLAGS) $(LIB) build/ch5md.o

lib: $(LIB)

clean:
	rm -f build/*.o $(LIB)

