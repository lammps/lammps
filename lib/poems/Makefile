# master -*- makefile -*- for linalg compilation

SHELL=/bin/sh

TARGETS= fermi 32-bit
CLEAN = $(patsubst %,clean-%,$(TARGETS))

default:
	@echo available targets: all $(TARGETS)

.PHONY: all default
.SUFFIXES:

all:
	$(MAKE) $(MFLAGS) $(TARGETS)

clean-all:
	$(MAKE) $(MFLAGS) $(CLEAN)

clean-%:
	$(MAKE) $(MFLAGS) -f $(patsubst clean-%,Makefile.%,$@) clean

%:
	$(MAKE) $(MFLAGS) -f Makefile.$@

