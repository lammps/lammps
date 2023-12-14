# Makefile for LAMMPS documentation

SHELL          = /bin/bash
HAS_BASH       = YES
ifeq (,$(wildcard $(SHELL)))
OSHELL         := $(SHELL)
override SHELL = /bin/sh
HAS_BASH       = NO
endif
BUILDDIR       = ${CURDIR}
RSTDIR         = $(BUILDDIR)/src
VENV           = $(BUILDDIR)/docenv
ANCHORCHECK    = $(VENV)/bin/rst_anchor_check
SPHINXCONFIG   = $(BUILDDIR)/utils/sphinx-config
MATHJAX        = $(SPHINXCONFIG)/_static/mathjax
MATHJAXTAG     = 3.2.2

PYTHON         = $(word 3,$(shell type python3))
DOXYGEN        = $(word 3,$(shell type doxygen))
HAS_PYTHON3    = NO
HAS_DOXYGEN    = NO
HAS_PDFLATEX   = NO

ifeq ($(shell type python3 >/dev/null 2>&1; echo $$?), 0)
HAS_PYTHON3    = YES
endif

ifeq ($(shell type doxygen >/dev/null 2>&1; echo $$?), 0)
HAS_DOXYGEN    = YES
endif

ifeq ($(shell type pdflatex >/dev/null 2>&1; echo $$?), 0)
ifeq ($(shell type latexmk >/dev/null 2>&1; echo $$?), 0)
HAS_PDFLATEX = YES
endif
endif

# override settings for PIP commands
# PIP_OPTIONS = --cert /etc/pki/ca-trust/extracted/openssl/ca-bundle.trust.crt --proxy http://proxy.mydomain.org

# temporarily disable caching so that the hack for the sphinx-tabs extensions to get proper non-html output works
SPHINXEXTRA = -j $(shell $(PYTHON) -c 'import multiprocessing;print(multiprocessing.cpu_count())')

# grab list of sources from doxygen config file.
# we only want to use explicitly listed files.
DOXYFILES      = $(shell sed -n -e 's/\#.*$$//' -e '/^ *INPUT \+=/,/^[A-Z_]\+ \+=/p' doxygen/Doxyfile.in | sed -e 's/@LAMMPS_SOURCE_DIR@/..\/src/g' -e 's/\\//g' -e 's/ \+/ /' -e 's/[A-Z_]\+ \+= *\(YES\|NO\|\)//') 

.PHONY: help clean-all clean clean-spelling epub mobi html pdf spelling anchor_check style_check char_check role_check xmlgen fasthtml

# ------------------------------------------

help:
	@if [ "$(HAS_BASH)" == "NO" ] ; then echo "bash was not found at $(OSHELL)! Please use: $(MAKE) SHELL=/path/to/bash" 1>&2; exit 1; fi
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html          create HTML pages in html dir"
	@echo "  pdf           create Manual.pdf in this dir"
	@echo "  epub          create ePUB format manual for e-book readers"
	@echo "  mobi          convert ePUB to MOBI format manual for e-book readers (e.g. Kindle)"
	@echo "                      (requires ebook-convert tool from calibre)"
	@echo "  fasthtml      approximate HTML page creation in fasthtml dir (for development)"
	@echo "  clean         remove all intermediate files"
	@echo "  clean-all     reset the entire build environment"
	@echo "  anchor_check  scan for duplicate anchor labels"
	@echo "  style_check   check for complete and consistent style lists"
	@echo "  package_check check for complete and consistent package lists"
	@echo "  role_check    check for misformatted role keywords"
	@echo "  spelling      spell-check the manual"

# ------------------------------------------

clean-all: clean
	rm -rf $(BUILDDIR)/docenv $(MATHJAX) $(BUILDDIR)/LAMMPS.mobi $(BUILDDIR)/LAMMPS.epub $(BUILDDIR)/Manual.pdf

clean: clean-spelling
	rm -rf $(BUILDDIR)/html $(BUILDDIR)/epub $(BUILDDIR)/latex $(BUILDDIR)/doctrees $(BUILDDIR)/doxygen/xml $(BUILDDIR)/doxygen-warn.log $(BUILDDIR)/doxygen/Doxyfile $(SPHINXCONFIG)/conf.py $(BUILDDIR)/fasthtml

clean-spelling:
	rm -rf $(BUILDDIR)/spelling

$(SPHINXCONFIG)/conf.py: $(SPHINXCONFIG)/conf.py.in
	sed -e 's,@DOXYGEN_XML_DIR@,$(BUILDDIR)/doxygen/xml,g'   \
	    -e 's,@LAMMPS_SOURCE_DIR@,$(BUILDDIR)/../src,g'    \
	    -e 's,@LAMMPS_PYTHON_DIR@,$(BUILDDIR)/../python,g' \
	    -e 's,@LAMMPS_DOC_DIR@,$(BUILDDIR),g' $< > $@

html: xmlgen $(VENV) $(SPHINXCONFIG)/conf.py $(ANCHORCHECK) $(MATHJAX)
	@if [ "$(HAS_BASH)" == "NO" ] ; then echo "bash was not found at $(OSHELL)! Please use: $(MAKE) SHELL=/path/to/bash" 1>&2; exit 1; fi
	@$(MAKE) $(MFLAGS) -C graphviz all
	@(\
		. $(VENV)/bin/activate ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		sphinx-build -E $(SPHINXEXTRA) -b html -c $(SPHINXCONFIG) -d $(BUILDDIR)/doctrees $(RSTDIR) html ;\
		touch $(RSTDIR)/Fortran.rst ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		sphinx-build $(SPHINXEXTRA) -b html -c $(SPHINXCONFIG) -d $(BUILDDIR)/doctrees $(RSTDIR) html ;\
		ln -sf Manual.html html/index.html;\
		rm -f $(BUILDDIR)/doxygen/xml/run.stamp;\
		echo "############################################" ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		rst_anchor_check src/*.rst ;\
		$(PYTHON) $(BUILDDIR)/utils/check-packages.py -s ../src -d src ;\
		env LC_ALL=C grep -n '[^ -~]' $(RSTDIR)/*.rst ;\
		env LC_ALL=C grep -n ' :[a-z]\+`' $(RSTDIR)/*.rst ;\
		env LC_ALL=C grep -n ' `[^`]\+<[a-z][^`]\+`[^_]' $(RSTDIR)/*.rst ;\
		env LC_ALL=C grep -n ':\(ref\|doc\):[^`]' $(RSTDIR)/*.rst ;\
		$(PYTHON) $(BUILDDIR)/utils/check-styles.py -s ../src -d src ;\
		echo "############################################" ;\
		deactivate ;\
	)
	@rm -rf html/_sources
	@rm -rf html/PDF
	@rm -rf html/USER
	@rm -rf html/JPG
	@cp -r src/PDF html/PDF
	@rm -rf html/PDF/.[sg]*
	@echo "Build finished. The HTML pages are in doc/html."

fasthtml: xmlgen $(VENV) $(SPHINXCONFIG)/conf.py $(ANCHORCHECK) $(MATHJAX)
	@if [ "$(HAS_BASH)" == "NO" ] ; then echo "bash was not found at $(OSHELL)! Please use: $(MAKE) SHELL=/path/to/bash" 1>&2; exit 1; fi
	@$(MAKE) $(MFLAGS) -C graphviz all
	@mkdir -p fasthtml
	@(\
		. $(VENV)/bin/activate ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		sphinx-build $(SPHINXEXTRA) -b html -c $(SPHINXCONFIG) -d $(BUILDDIR)/fasthtml/doctrees $(RSTDIR) fasthtml ;\
		touch $(RSTDIR)/Fortran.rst ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		sphinx-build $(SPHINXEXTRA) -b html -c $(SPHINXCONFIG) -d $(BUILDDIR)/fasthtml/doctrees $(RSTDIR) fasthtml ;\
		deactivate ;\
	)
	@rm -rf fasthtml/_sources
	@rm -rf fasthtml/PDF
	@rm -rf fasthtml/USER
	@rm -rf fasthtml/JPG
	@cp -r src/PDF fasthtml/PDF
	@rm -rf fasthtml/PDF/.[sg]*
	@echo "Fast HTML build finished. The HTML pages are in doc/fasthtml."

spelling: xmlgen $(SPHINXCONFIG)/conf.py $(VENV) $(SPHINXCONFIG)/false_positives.txt
	@if [ "$(HAS_BASH)" == "NO" ] ; then echo "bash was not found at $(OSHELL)! Please use: $(MAKE) SHELL=/path/to/bash" 1>&2; exit 1; fi
	@(\
		. $(VENV)/bin/activate ; \
		cp $(SPHINXCONFIG)/false_positives.txt $(RSTDIR)/;  env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		sphinx-build -b spelling -c $(SPHINXCONFIG) -d $(BUILDDIR)/doctrees $(RSTDIR) spelling ;\
		rm -f $(BUILDDIR)/doxygen/xml/run.stamp;\
		deactivate ;\
	)
	@echo "Spell check finished."

epub: xmlgen $(VENV) $(SPHINXCONFIG)/conf.py $(ANCHORCHECK)
	@if [ "$(HAS_BASH)" == "NO" ] ; then echo "bash was not found at $(OSHELL)! Please use: $(MAKE) SHELL=/path/to/bash" 1>&2; exit 1; fi
	@$(MAKE) $(MFLAGS) -C graphviz all
	@mkdir -p epub/JPG
	@rm -f LAMMPS.epub
	@cp src/JPG/*.* epub/JPG
	@(\
		. $(VENV)/bin/activate ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		sphinx-build -E $(SPHINXEXTRA) -b epub -c $(SPHINXCONFIG) -d $(BUILDDIR)/doctrees $(RSTDIR) epub ;\
		touch $(RSTDIR)/Fortran.rst ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		sphinx-build $(SPHINXEXTRA) -b epub -c $(SPHINXCONFIG) -d $(BUILDDIR)/doctrees $(RSTDIR) epub ;\
		rm -f $(BUILDDIR)/doxygen/xml/run.stamp;\
		deactivate ;\
	)
	@mv  epub/LAMMPS.epub .
	@rm -rf epub
	@echo "Build finished. The ePUB manual file is created."

mobi: epub
	@rm -f LAMMPS.mobi
	@ebook-convert LAMMPS.epub LAMMPS.mobi
	@echo "Conversion finished. The MOBI manual file is created."

pdf: xmlgen $(VENV) $(SPHINXCONFIG)/conf.py $(ANCHORCHECK)
	@if [ "$(HAS_BASH)" == "NO" ] ; then echo "bash was not found at $(OSHELL)! Please use: $(MAKE) SHELL=/path/to/bash" 1>&2; exit 1; fi
	@$(MAKE) $(MFLAGS) -C graphviz all
	@if [ "$(HAS_PDFLATEX)" == "NO" ] ; then echo "PDFLaTeX or latexmk were not found! Please check README for further instructions" 1>&2; exit 1; fi
	@(\
		. $(VENV)/bin/activate ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		sphinx-build -E $(SPHINXEXTRA) -b latex -c $(SPHINXCONFIG) -d $(BUILDDIR)/doctrees $(RSTDIR) latex ;\
		touch $(RSTDIR)/Fortran.rst ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		sphinx-build  $(SPHINXEXTRA) -b latex -c $(SPHINXCONFIG) -d $(BUILDDIR)/doctrees $(RSTDIR) latex ;\
		rm -f $(BUILDDIR)/doxygen/xml/run.stamp;\
		echo "############################################" ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		rst_anchor_check src/*.rst ;\
		$(PYTHON) utils/check-packages.py -s ../src -d src ;\
		env LC_ALL=C grep -n '[^ -~]' $(RSTDIR)/*.rst ;\
		env LC_ALL=C grep -n ' :[a-z]\+`' $(RSTDIR)/*.rst ;\
		env LC_ALL=C grep -n ' `[^`]\+<[a-z][^`]\+`[^_]' $(RSTDIR)/*.rst ;\
		env LC_ALL=C grep -n ':\(ref\|doc\):[^`]' $(RSTDIR)/*.rst ;\
		$(PYTHON) utils/check-styles.py -s ../src -d src ;\
		echo "############################################" ;\
		deactivate ;\
	)
	@cd latex && \
		sed 's/\\begin{equation}//g' LAMMPS.tex > tmp.tex && \
		mv tmp.tex LAMMPS.tex && \
		sed 's/\\end{equation}//g' LAMMPS.tex > tmp.tex && \
		mv tmp.tex LAMMPS.tex && \
		sed 's/\\contentsname}{.*}}/\\contentsname}{LAMMPS Documentation}}/g' LAMMPS.tex > tmp.tex && \
		mv tmp.tex LAMMPS.tex && \
		$(MAKE) $(MFLAGS) && \
		mv LAMMPS.pdf ../Manual.pdf && \
		cd ../;
	@rm -rf latex/_sources
	@rm -rf latex/PDF
	@rm -rf latex/USER
	@cp -r src/PDF latex/PDF
	@rm -rf latex/PDF/.[sg]*
	@echo "Build finished. Manual.pdf is in this directory."

anchor_check : $(ANCHORCHECK)
	@(\
		. $(VENV)/bin/activate ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		rst_anchor_check src/*.rst ;\
		deactivate ;\
	)

style_check : $(VENV)
	@(\
		. $(VENV)/bin/activate ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		$(PYTHON) utils/check-styles.py -s ../src -d src ;\
		deactivate ;\
	)

package_check : $(VENV)
	@(\
		. $(VENV)/bin/activate ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		$(PYTHON) utils/check-packages.py -s ../src -d src ;\
		deactivate ;\
	)

char_check :
	@( env LC_ALL=C grep -n '[^ -~]' $(RSTDIR)/*.rst && exit 1 || : )

role_check :
	@( env LC_ALL=C grep -n ' :[a-z]\+`' $(RSTDIR)/*.rst && exit 1 || : )
	@( env LC_ALL=C grep -n ' `[^`]\+<[a-z][^`]\+`[^_]' $(RSTDIR)/*.rst && exit 1 || : )
	@( env LC_ALL=C grep -n ':\(ref\|doc\):[^`]' $(RSTDIR)/*.rst && exit 1 || : )

link_check : $(VENV) html
	@(\
		. $(VENV)/bin/activate ; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		linkchecker -F html --check-extern html/Manual.html ;\
		deactivate ;\
	)

xmlgen : doxygen/xml/index.xml

doxygen/Doxyfile: doxygen/Doxyfile.in
	sed -e 's/@LAMMPS_SOURCE_DIR@/..\/..\/src/g' $< > $@

doxygen/xml/index.xml : $(VENV) doxygen/Doxyfile $(DOXYFILES)
	@(cd doxygen; $(DOXYGEN) Doxyfile && touch xml/run.stamp)
# ------------------------------------------

$(VENV):
	@if [ "$(HAS_BASH)" == "NO" ] ; then echo "bash was not found at $(OSHELL)! Please use: $(MAKE) SHELL=/path/to/bash" 1>&2; exit 1; fi
	@if [ "$(HAS_PYTHON3)" == "NO" ] ; then echo "python3 was not found! Please see README for further instructions" 1>&2; exit 1; fi
	@if [ "$(HAS_DOXYGEN)" == "NO" ] ; then echo "doxygen was not found! Please see README for further instructions" 1>&2; exit 1; fi
	@( \
		$(PYTHON) -m venv $(VENV); \
		. $(VENV)/bin/activate; \
		pip $(PIP_OPTIONS) install --upgrade pip; \
		pip $(PIP_OPTIONS) install --upgrade wheel; \
		pip $(PIP_OPTIONS) install -r $(BUILDDIR)/utils/requirements.txt; \
		deactivate;\
	)

$(MATHJAX):
	@git clone -b $(MATHJAXTAG) -c advice.detachedHead=0 --depth 1 https://github.com/mathjax/MathJax.git $@

$(ANCHORCHECK): $(VENV)
	@( \
		. $(VENV)/bin/activate; env PYTHONWARNINGS= PYTHONDONTWRITEBYTECODE=1 \
		pip $(PIP_OPTIONS) install -e utils/converters;\
		deactivate;\
	)
