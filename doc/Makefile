# Makefile for LAMMPS documentation

SHELL          = /bin/bash
BUILDDIR       = ${CURDIR}
RSTDIR         = $(BUILDDIR)/src
VENV           = $(BUILDDIR)/docenv
MATHJAX        = $(BUILDDIR)/mathjax
TXT2RST        = $(VENV)/bin/txt2rst
ANCHORCHECK    = $(VENV)/bin/rst_anchor_check

PYTHON         = $(shell which python3)
DOXYGEN        = $(shell which doxygen)
VIRTUALENV     = virtualenv
HAS_PYTHON3    = NO
HAS_VIRTUALENV = NO
HAS_DOXYGEN    = NO

ifeq ($(shell which python3 >/dev/null 2>&1; echo $$?), 0)
HAS_PYTHON3    = YES
endif

ifeq ($(shell which doxygen >/dev/null 2>&1; echo $$?), 0)
HAS_DOXYGEN    = YES
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

# grab list of sources from doxygen config file.
# we only want to use explicitly listed files.
DOXYFILES      = $(shell sed -n -e 's/\#.*$$//' -e '/^ *INPUT \+=/,/^[A-Z_]\+ \+=/p' doxygen/Doxyfile | sed -e 's/\.\.\/src/src/g' -e 's/\\//g' -e 's/ \+/ /' -e 's/[A-Z_]\+ \+= *\(YES\|NO\|\)//') 

.PHONY: help clean-all clean clean-spelling epub mobi rst html pdf spelling anchor_check style_check xmlgen

# ------------------------------------------

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html          create HTML doc pages in html dir"
	@echo "  pdf           create Developer.pdf and Manual.pdf in this dir"
	@echo "  fetch         fetch HTML and PDF files from LAMMPS web site"
	@echo "  epub          create ePUB format manual for e-book readers"
	@echo "  mobi          convert ePUB to MOBI format manual for e-book readers (e.g. Kindle)"
	@echo "                      (requires ebook-convert tool from calibre)"
	@echo "  clean         remove all intermediate RST files"
	@echo "  clean-all     reset the entire build environment"
	@echo "  anchor_check  scan for duplicate anchor labels"
	@echo "  style_check   check for complete and consistent style lists"
	@echo "  package_check check for complete and consistent package lists"
	@echo "  spelling      spell-check the manual"

# ------------------------------------------

clean-all: clean
	rm -rf $(BUILDDIR)/docenv $(BUILDDIR)/mathjax LAMMPS.mobi LAMMPS.epub Manual.pdf Developer.pdf

clean: clean-spelling
	rm -rf html $(BUILDDIR)/epub $(BUILDDIR)/latex $(BUILDDIR)/doctrees $(BUILDDIR)/doxygen/xml $(BUILDDIR)/doxygen-warn.log

clean-spelling:
	rm -rf spelling

html: xmlgen $(ANCHORCHECK) $(MATHJAX)
	@(\
		. $(VENV)/bin/activate ;\
		sphinx-build $(SPHINXEXTRA) -b html -c utils/sphinx-config -d $(BUILDDIR)/doctrees $(RSTDIR) html ;\
		echo "############################################" ;\
		rst_anchor_check src/*.rst ;\
		python utils/check-packages.py -s ../src -d src ;\
		env LC_ALL=C grep -n '[^ -~]' $(RSTDIR)/*.rst ;\
		python utils/check-styles.py -s ../src -d src ;\
		echo "############################################" ;\
		deactivate ;\
	)
	@rm -rf html/_sources
	@rm -rf html/PDF
	@rm -rf html/USER
	@rm -rf html/JPG
	@cp -r src/PDF html/PDF
	@mkdir -p html/JPG
	@cp `grep -A2 '\.\. .*image::' src/*.rst | grep ':target:' | sed -e 's,.*:target: JPG/,src/JPG/,' | sort | uniq` html/JPG/
	@rm -rf html/PDF/.[sg]*
	@mkdir -p html/_static/mathjax
	@cp -r $(MATHJAX)/es5 html/_static/mathjax/
	@echo "Build finished. The HTML pages are in doc/html."

spelling: $(VENV) utils/sphinx-config/false_positives.txt
	@(\
		. $(VENV)/bin/activate ;\
		pip install sphinxcontrib-spelling ;\
		cp utils/sphinx-config/false_positives.txt $(RSTDIR)/ ;\
		sphinx-build -b spelling -c utils/sphinx-config -d $(BUILDDIR)/doctrees $(RSTDIR) spelling ;\
		deactivate ;\
	)
	@echo "Spell check finished."

epub: xmlgen $(VENV)
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

pdf: xmlgen $(ANCHORCHECK)
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
		python utils/check-packages.py -s ../src -d src ;\
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
	@rm -rf latex/PDF/.[sg]*
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

style_check : $(VENV)
	@(\
		. $(VENV)/bin/activate ;\
		python utils/check-styles.py -s ../src -d src ;\
		deactivate ;\
	)

package_check : $(VENV)
	@(\
		. $(VENV)/bin/activate ;\
		python utils/check-packages.py -s ../src -d src ;\
		deactivate ;\
	)

xmlgen : doxygen/xml/index.xml

doxygen/xml/index.xml : $(VENV) doxygen/Doxyfile $(DOXYFILES)
	@(cd doxygen; doxygen Doxyfile)
	@touch src/pg_*.rst

# ------------------------------------------

$(VENV):
	@if [ "$(HAS_PYTHON3)" == "NO" ] ; then echo "python3 was not found! Please see README for further instructions" 1>&2; exit 1; fi
	@if [ "$(HAS_DOXYGEN)" == "NO" ] ; then echo "doxygen was not found! Please see README for further instructions" 1>&2; exit 1; fi
	@if [ "$(HAS_VIRTUALENV)" == "NO" ] ; then echo "virtualenv was not found! Please see README for further instructions" 1>&2; exit 1; fi
	@( \
		$(VIRTUALENV) -p $(PYTHON) $(VENV); \
		. $(VENV)/bin/activate; \
		pip install Sphinx; \
		pip install breathe; \
		deactivate;\
	)

$(MATHJAX):
	@git clone --depth 1 https://github.com/mathjax/MathJax.git mathjax

$(TXT2RST) $(ANCHORCHECK): $(VENV)
	@( \
		. $(VENV)/bin/activate; \
		(cd utils/converters;\
		python setup.py develop);\
		deactivate;\
	)
