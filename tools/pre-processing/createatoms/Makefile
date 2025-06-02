# Makefile for the chain tools

#
# Targets
#

all:
	$(MAKE) createAtoms

createAtoms:	createAtoms.f
	gfortran createAtoms.f -o createAtoms

clean:
	rm createAtoms

#
# Rules
#

.f.o:	
	gfortran -O -w -c $<
