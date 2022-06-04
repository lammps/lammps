# Makefile for serial tools

#
# Targets
#

all:
	$(MAKE) binary2txt chain micelle2d stl_bin2txt

binary2txt:	binary2txt.o
	g++ -g binary2txt.o -o binary2txt

chain:	chain.o
	ifort chain.o -o chain

micelle2d:	micelle2d.o
	ifort micelle2d.o -o micelle2d

thermo_extract:	thermo_extract.o
	gcc -g thermo_extract.o -o thermo_extract

stl_bin2txt:	stl_bin2txt.o
	g++ -g stl_bin2txt.o -o stl_bin2txt

clean:
	rm binary2txt chain micelle2d
	rm thermo_extract
	rm *.o

#
# Rules
#

.cpp.o:
	g++ -g -c $<

.c.o:
	gcc -g -c $<

.f.o:	
	ifort -O -w -c $<
