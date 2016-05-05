# Makefile for LAMMPS documentation
SHA1          = $(shell echo $USER-$PWD | sha1sum | cut -f1 -d" ")
BUILDDIR      = /tmp/lammps-docs-$(SHA1)
RSTDIR        = $(BUILDDIR)/rst
VENV          = $(BUILDDIR)/docenv
TXT2RST       = $(VENV)/bin/txt2rst

PYTHON        = $(shell which python3)

ifeq ($(shell which python3 >/dev/null 2>&1; echo $$?), 1)
$(error Python3 was not found! Please check README.md for further instructions)
endif

ifeq ($(shell which virtualenv >/dev/null 2>&1; echo $$?), 1)
$(error virtualenv was not found! Please check README.md for further instructions)
endif

SOURCES=$(wildcard src/*.txt)
OBJECTS=$(SOURCES:src/%.txt=$(RSTDIR)/%.rst)

.PHONY: help clean-all clean html pdf venv

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make HTML version of documentation using Sphinx"
	@echo "  pdf        to make Manual.pdf"
	@echo "  clean      to remove all generated RST files"
	@echo "  clean-all  to reset the entire build environment"

clean-all:
	rm -rf $(BUILDDIR)/*

clean:
	rm -rf $(RSTDIR)

html: $(OBJECTS)
	@(\
		source $(VENV)/bin/activate ;\
		cp -r src/* $(RSTDIR)/ ;\
		sphinx-build -j 8 -b html -c utils/sphinx-config -d $(BUILDDIR)/doctrees $(RSTDIR) html ;\
		deactivate ;\
	)	
	@echo "Build finished. The HTML pages are in doc/html."

pdf: html
	htmldoc --title --toctitle "Table of Contents" --tocfooter ..i --toclevels 4 --header ... --footer ..1 --size letter --linkstyle plain --linkcolor blue -f Manual.pdf html/Manual.html html/Section_intro.html html/Section_start.html html/Section_commands.html html/Section_packages.html html/Section_accelerate.html html/Section_howto.html html/Section_example.html html/Section_perf.html html/Section_tools.html html/Section_modify.html html/Section_python.html html/Section_errors.html html/Section_history.html html/[a-z]*.html

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
