# Makefile for LAMMPS documentation

SHELL         = /bin/bash
BUILDDIR      = ${CURDIR}
RSTDIR        = $(BUILDDIR)/src
VENV          = $(BUILDDIR)/docenv
TXT2RST       = $(VENV)/bin/txt2rst
ANCHORCHECK   = $(VENV)/bin/rst_anchor_check

PYTHON        = $(shell which python3)
VIRTUALENV     = virtualenv
HAS_PYTHON3    = NO
HAS_VIRTUALENV = NO

ifeq ($(shell which python3 >/dev/null 2>&1; echo $$?), 0)
HAS_PYTHON3 = YES
endif

ifeq ($(shell which virtualenv-3 >/dev/null 2>&1; echo $$?), 0)
VIRTUALENV     = virtualenv-3
HAS_VIRTUALENV = YES
endif

ifeq ($(shell which virtualenv >/dev/null 2>&1; echo $$?), 0)
VIRTUALENV     = virtualenv
HAS_VIRTUALENV = YES
endif

SPHINXEXTRA = -j $(shell $(PYTHON) -c 'import multiprocessing;print(multiprocessing.cpu_count())')

.PHONY: help clean-all clean epub mobi rst html pdf venv spelling anchor_check style_check

# ------------------------------------------

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       create HTML doc pages in html dir"
	@echo "  pdf        create Developer.pdf and Manual.pdf in this dir"
	@echo "  fetch      fetch HTML and PDF files from LAMMPS web site"
	@echo "  epub       create ePUB format manual for e-book readers"
	@echo "  mobi       convert ePUB to MOBI format manual for e-book readers (e.g. Kindle)"
	@echo "                      (requires ebook-convert tool from calibre)"
	@echo "  clean      remove all intermediate RST files"
	@echo "  clean-all  reset the entire build environment"
	@echo "  anchor_check  scan for duplicate anchor labels"
	@echo "  style_check   check for complete and consistent style lists"
	@echo "  spelling   spell-check the manual"

# ------------------------------------------

clean-all: clean
	rm -rf $(BUILDDIR)/docenv $(BUILDDIR)/doctrees

clean:
	rm -rf html epub latex
	rm -rf spelling

clean-spelling:
	rm -rf spelling

rst: clean $(ANCHORCHECK)

html: $(ANCHORCHECK)
	@(\
		. $(VENV)/bin/activate ;\
		sphinx-build $(SPHINXEXTRA) -b html -c utils/sphinx-config -d $(BUILDDIR)/doctrees $(RSTDIR) html ;\
		echo "############################################" ;\
		rst_anchor_check src/*.rst ;\
		env LC_ALL=C grep -n '[^ -~]' $(RSTDIR)/*.rst ;\
		python utils/check-styles.py -s ../src -d src ;\
		echo "############################################" ;\
		deactivate ;\
	)
	-rm html/searchindex.js
	@rm -rf html/_sources
	@rm -rf html/PDF
	@rm -rf html/USER
	@rm -rf html/JPG
	@cp -r src/PDF html/PDF
	@cp -r src/USER html/USER
	@mkdir -p html/JPG
	@cp `grep -A2 '\.\. image::' src/*.rst | grep ':target:' | sed -e 's,.*:target: JPG/,src/JPG/,' | sort | uniq` html/JPG/
	@rm -rf html/PDF/.[sg]*
	@rm -rf html/USER/.[sg]*
	@rm -rf html/USER/*/.[sg]*
	@rm -rf html/USER/*/*.[sg]*
	@echo "Build finished. The HTML pages are in doc/html."

spelling: utils/sphinx-config/false_positives.txt
	@(\
		. $(VENV)/bin/activate ;\
		pip install sphinxcontrib-spelling ;\
		cp utils/sphinx-config/false_positives.txt $(RSTDIR)/ ;\
		sphinx-build -b spelling -c utils/sphinx-config -d $(BUILDDIR)/doctrees $(RSTDIR) spelling ;\
		deactivate ;\
	)
	@echo "Spell check finished."

epub:
	@mkdir -p epub/JPG
	@rm -f LAMMPS.epub
	@cp src/JPG/lammps-logo.png epub/
	@cp src/JPG/*.* epub/JPG
	@(\
		. $(VENV)/bin/activate ;\
		sphinx-build $(SPHINXEXTRA) -b epub -c utils/sphinx-config -d $(BUILDDIR)/doctrees $(RSTDIR) epub ;\
		deactivate ;\
	)
	@mv  epub/LAMMPS.epub .
	@rm -rf epub
	@echo "Build finished. The ePUB manual file is created."

mobi: epub
	@rm -f LAMMPS.mobi
	@ebook-convert LAMMPS.epub LAMMPS.mobi
	@echo "Conversion finished. The MOBI manual file is created."

pdf: $(ANCHORCHECK)
	@(\
		cd src/Developer; \
		pdflatex developer; \
		pdflatex developer; \
		mv developer.pdf ../../Developer.pdf; \
		cd ../../; \
	)
	@(\
		. $(VENV)/bin/activate ;\
		sphinx-build $(SPHINXEXTRA) -b latex -c utils/sphinx-config -d $(BUILDDIR)/doctrees $(RSTDIR) latex ;\
		echo "############################################" ;\
		rst_anchor_check src/*.rst ;\
		env LC_ALL=C grep -n '[^ -~]' $(RSTDIR)/*.rst ;\
		python utils/check-styles.py -s ../src -d src ;\
		echo "############################################" ;\
		deactivate ;\
	)
	@cd latex && \
		sed 's/latexmk -pdf -dvi- -ps-/pdflatex/g' Makefile > temp && \
		mv temp Makefile && \
		sed 's/\\begin{equation}//g' LAMMPS.tex > tmp.tex && \
		mv tmp.tex LAMMPS.tex && \
		sed 's/\\end{equation}//g' LAMMPS.tex > tmp.tex && \
		mv tmp.tex LAMMPS.tex && \
		make && \
		make && \
		make && \
		mv LAMMPS.pdf ../Manual.pdf && \
		cd ../;
	@rm -rf latex/_sources
	@rm -rf latex/PDF
	@rm -rf latex/USER
	@cp -r src/PDF latex/PDF
	@cp -r src/USER latex/USER
	@rm -rf latex/PDF/.[sg]*
	@rm -rf latex/USER/.[sg]*
	@rm -rf latex/USER/*/.[sg]*
	@rm -rf latex/USER/*/*.[sg]*
	@echo "Build finished. Manual.pdf and Developer.pdf are in this directory."

fetch:
	@rm -rf html_www Manual_www.pdf Developer_www.pdf
	@curl -s -o Manual_www.pdf http://lammps.sandia.gov/doc/Manual.pdf
	@curl -s -o Developer_www.pdf http://lammps.sandia.gov/doc/Developer.pdf
	@curl -s -o lammps-doc.tar.gz http://lammps.sandia.gov/tars/lammps-doc.tar.gz
	@tar xzf lammps-doc.tar.gz
	@rm -f lammps-doc.tar.gz

anchor_check : $(ANCHORCHECK)
	@(\
		. $(VENV)/bin/activate ;\
		rst_anchor_check src/*.rst ;\
		deactivate ;\
	)

style_check :
	@(\
		. $(VENV)/bin/activate ;\
		python utils/check-styles.py -s ../src -d src ;\
		deactivate ;\
	)

# ------------------------------------------

$(VENV):
	@if [ "$(HAS_PYTHON3)" == "NO" ] ; then echo "Python3 was not found! Please check README.md for further instructions" 1>&2; exit 1; fi
	@if [ "$(HAS_VIRTUALENV)" == "NO" ] ; then echo "virtualenv was not found! Please check README.md for further instructions" 1>&2; exit 1; fi
	@( \
		$(VIRTUALENV) -p $(PYTHON) $(VENV); \
		. $(VENV)/bin/activate; \
		pip install Sphinx; \
		deactivate;\
	)

$(TXT2RST) $(ANCHORCHECK): $(VENV)
	@( \
		. $(VENV)/bin/activate; \
		(cd utils/converters;\
		python setup.py develop);\
		deactivate;\
	)
