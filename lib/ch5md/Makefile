EXTRAMAKE=Makefile.lammps.empty

CC=h5cc
# -DH5_NO_DEPRECATED_SYMBOLS is required to force using the v1.8 API
# on HDF5 installations configured to default to using the v1.6 API.
CFLAGS=-D_DEFAULT_SOURCE -O2 -DH5_NO_DEPRECATED_SYMBOLS 
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

