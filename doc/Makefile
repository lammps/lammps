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
	@echo "  txt2html   to build txt2html tool"
	@echo "  clean      to remove all generated RST files"
	@echo "  clean-all  to reset the entire build environment"

clean-all:
	rm -rf $(BUILDDIR)/* utils/txt2html/txt2html.exe

clean:
	rm -rf $(RSTDIR)

txt2html: utils/txt2html/txt2html.exe

html: $(OBJECTS)
	@(\
		. $(VENV)/bin/activate ;\
		cp -r src/* $(RSTDIR)/ ;\
		sphinx-build -j 8 -b html -c utils/sphinx-config -d $(BUILDDIR)/doctrees $(RSTDIR) html ;\
		deactivate ;\
	)	
	@echo "Build finished. The HTML pages are in doc/html."

pdf: utils/txt2html/txt2html.exe
	@(\
		cd src; \
		../utils/txt2html/txt2html.exe -b *.txt; \
		htmldoc --batch ../lammps.book;          \
		for s in `echo *.txt | sed -e 's,\.txt,\.html,g'` ; \
			do grep -q $$s ../lammps.book || \
			echo doc file $$s missing in lammps.book; done; \
		rm *.html; \
	)

utils/txt2html/txt2html.exe: utils/txt2html/txt2html.cpp
	g++ -O -Wall -o $@ $<

$(RSTDIR)/%.rst : src/%.txt $(TXT2RST)
	@(\
		mkdir -p $(RSTDIR) ; \
		. $(VENV)/bin/activate ;\
		txt2rst $< > $@ ;\
		deactivate ;\
	)

$(VENV):
	@( \
		virtualenv -p $(PYTHON) $(VENV); \
		. $(VENV)/bin/activate; \
		pip install Sphinx; \
		pip install sphinxcontrib-images; \
		deactivate;\
	)

$(TXT2RST): $(VENV)
	@( \
		. $(VENV)/bin/activate; \
		(cd utils/converters;\
		python setup.py develop);\
		deactivate;\
	)
