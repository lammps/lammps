# master -*- makefile -*- for libqmmm compilation

SHELL=/bin/sh

TARGETS= gfortran
CLEAN = $(patsubst %,clean-%,$(TARGETS))

default:
	@echo available targets: $(TARGETS)

.PHONY: default
.SUFFIXES:

all:
	$(MAKE) $(MFLAGS) $(TARGETS)

clean clean-all:
	$(MAKE) $(MFLAGS) $(CLEAN)

clean-%:
	$(MAKE) $(MFLAGS) -f $(patsubst clean-%,Makefile.%,$@) clean

%:
	$(MAKE) $(MFLAGS) -f Makefile.$@

