EXTRAMAKE=Makefile.lammps.empty

CC=h5cc
HDF5_PATH=/usr
INC=-I include
AR=ar
ARFLAGS=rc
LIB=libch5md.a

all: lib

ch5md.o: src/ch5md.c
	$(CC) $(INC) -c $< -o build/$@

lib: ch5md.o
	$(AR) $(ARFLAGS) $(LIB) build/ch5md.o

clean:
	rm -f build/*.o $(LIB)

ch5md_SYSINC=$(HDF5_PATH)/include
ch5md_SYSLIB=$(HDF5_PATH)/lib
