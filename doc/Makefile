# Makefile for LAMMPS documentation
BUILDDIR      = /tmp/lammps-docs
RSTDIR        = /tmp/lammps-docs/rst
VENV          = $(BUILDDIR)/docenv
TXT2RST       = $(VENV)/bin/txt2rst

PYTHON        = $(shell which python3)

ifeq ($(shell which python3 >/dev/null 2>&1; echo $$?), 1)
$(error Python3 was not found)
endif

ifeq ($(shell which virtualenv >/dev/null 2>&1; echo $$?), 1)
$(error virtualenv was not found)
endif

SOURCES=$(wildcard src/*.txt)
OBJECTS=$(SOURCES:src/%.txt=$(RSTDIR)/%.rst)

.PHONY: help clean html pdf venv

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make HTML version of documentation using Sphinx"
	@echo "  pdf        to make Manual.pdf"

clean:
	rm -rf $(BUILDDIR)/*

html: $(OBJECTS)
	@(\
		source $(VENV)/bin/activate ;\
		cp -r src/* $(RSTDIR)/ ;\
		sphinx-build -j 8 -b html -c utils/sphinx-config -d $(BUILDDIR)/doctrees $(RSTDIR) html ;\
		deactivate ;\
	)	
	@echo "Build finished. The HTML pages are in doc/html."


$(RSTDIR)/%.rst : src/%.txt $(TXT2RST)
	@(\
		mkdir -p $(RSTDIR) ; \
		source $(VENV)/bin/activate ;\
		txt2rst $< > $@ ;\
		deactivate ;\
	)

$(VENV):
	@( \
		virtualenv -p $(PYTHON) $(VENV); \
		source $(VENV)/bin/activate; \
		pip install Sphinx; \
		pip install sphinxcontrib-images; \
		deactivate;\
	)

$(TXT2RST): $(VENV)
	@( \
		source $(VENV)/bin/activate; \
		pushd utils/converters;\
		python setup.py develop;\
		popd;\
		deactivate;\
	)
