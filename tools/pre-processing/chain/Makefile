# Makefile for the chain tools

#
# Targets
#

all:
	$(MAKE) chain

chain:	chain.f90
	gfortran chain.f90 -o chain

clean:
	rm chain
	rm *.mod

#
# Rules
#

.f.o:	
	gfortran -O -w -c $<
