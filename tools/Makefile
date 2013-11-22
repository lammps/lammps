# Makefile for serial tools

#
# Targets
#

all:
	$(MAKE) binary2txt restart2data chain micelle2d data2xmovie

binary2txt:	binary2txt.o
	g++ -g binary2txt.o -o binary2txt

chain:	chain.o
	ifort chain.o -o chain

micelle2d:	micelle2d.o
	ifort micelle2d.o -o micelle2d

data2xmovie:	data2xmovie.o
	g++ -g data2xmovie.o -o data2xmovie

thermo_extract:	thermo_extract.o
	gcc -g thermo_extract.o -o thermo_extract

clean:
	rm binary2txt restart2data chain micelle2d data2xmovie
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
