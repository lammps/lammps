EXTRAMAKE=Makefile.lammps.empty

CC=h5cc
HDF5_PATH=/usr
INC=-I include
AR=ar
ARFLAGS=rc
LIB=libch5md.a

all: lib

build:
	mkdir -p build
build/ch5md.o: src/ch5md.c | build
	$(CC) $(INC) -c $< -o $@

.PHONY: Makefile.lammps
Makefile.lammps:
	@echo "# Settings that the LAMMPS build will import when this package library is used" > Makefile.lammps
	@echo "ch5md_SYSINC = -I${HDF5_PATH}/include" >> Makefile.lammps
	@echo "ch5md_SYSLIB = -L${HDF5_PATH}/lib" >> Makefile.lammps

$(LIB): build/ch5md.o
	$(AR) $(ARFLAGS) $(LIB) build/ch5md.o

lib: $(LIB) Makefile.lammps

clean:
	rm -f build/*.o $(LIB) Makefile.lammps

