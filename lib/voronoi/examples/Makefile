# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Makefile rules
all: ex_basic ex_walls ex_custom ex_extra ex_degenerate ex_interface

ex_basic:
	$(MAKE) -C basic

ex_walls:
	$(MAKE) -C walls

ex_custom:
	$(MAKE) -C custom

ex_extra:
	$(MAKE) -C extra

ex_degenerate:
	$(MAKE) -C degenerate

ex_interface:
	$(MAKE) -C interface

clean:
	$(MAKE) -C basic clean
	$(MAKE) -C walls clean
	$(MAKE) -C custom clean
	$(MAKE) -C extra clean
	$(MAKE) -C degenerate clean
	$(MAKE) -C interface clean

.PHONY: all ex_basic ex_walls ex_custom ex_extra ex_degenerate clean
