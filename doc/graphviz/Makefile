# Makefile for generating images with graphviz
#
BUILDDIR   = ${CURDIR}/..
IMGDIR     = $(BUILDDIR)/src/JPG
IMGSRC     = $(wildcard *.dot)
IMGPNG     = $(IMGSRC:%.dot=$(IMGDIR)/%.png)

HAS_DOT        = NO
ifeq ($(shell which dot >/dev/null 2>&1; echo $$?), 0)
HAS_DOT        = YES
endif

all:    $(IMGPNG)

clean:
	rm -f $(IMGSVG) $(IMGPDF) $(IMGPNG) *~

ifeq ($(HAS_DOT),YES)
$(IMGDIR)/lammps-classes.png : lammps-classes.dot
	dot -Tpng -Kneato -o $@ $<

$(IMGDIR)/%.png: %.dot
	dot -Tpng -Kdot -o $@ $<
endif

ifeq ($(HAS_DOT),NO)
$(IMGDIR)/%.png: %.dot
	@echo '###################################################'
	@echo '# Need to install "graphviz" to regenerate graphs #'
	@echo '###################################################'
endif

