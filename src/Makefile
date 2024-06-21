# LAMMPS multiple-machine -*- Makefile -*-

SHELL = /bin/bash
PYTHON = python
DYN_LIB = -ldl

#.IGNORE:

# Definitions

ROOT =   lmp
EXE =    lmp_$@
ARLIB =  liblammps_$@.a
SHLIB =  liblammps_$@.so
ARLINK = liblammps.a
SHLINK = liblammps.so
TMPNAME= tmp_$@_name
LMPLINK= -L. -llammps_$@

OBJDIR =   Obj_$@
OBJSHDIR = Obj_shared_$@

SRC =   $(wildcard *.cpp)
INC =   $(filter-out lmpinstalledpkgs.h lmpgitversion.h,$(wildcard *.h))
OBJ =   $(SRC:.cpp=.o)

SRCLIB = $(filter-out main.cpp,$(SRC))
OBJLIB = $(filter-out main.o,$(OBJ))

# Command-line options for mode: static (default), shared, or print

mode = static
objdir = $(OBJDIR)

ifeq ($(mode),shared)
objdir = $(OBJSHDIR)
endif

# Package variables

# PACKAGE    = all packages
# PACKBASIC  = a few core packages
# PACKMOST   = most packages that do not require additional libs
# PACKLIB    = all packages that require an additional lib
#              should be PACKSYS + PACKINT + PACKEXT
# PACKSYS    = subset that reqiure a common system library
#              include LATBOLTZ b/c it requires a full MPI, not just STUBS
# PACKINT    = subset that require an internal (provided) library
# PACKEXT    = subset that require an external (downloaded) library

PACKAGE = \
	molecule \
	kspace \
	adios \
	amoeba \
	asphere \
	awpmd \
	bocs \
	body \
	bpm \
	brownian \
	cg-dna \
	cg-spica \
	class2 \
	colloid \
	colvars \
	compress \
	coreshell \
	diffraction \
	dipole \
	dpd-basic \
	dpd-meso \
	dpd-react \
	dpd-smooth \
	drude \
	eff \
	extra-command \
	extra-compute \
	extra-dump \
	extra-fix \
	extra-molecule \
	extra-pair \
	fep \
	gpu \
	granular \
	h5md \
	intel \
	interlayer \
	kim \
	kokkos \
	latboltz \
	lepton \
	machdyn \
	manifold \
	manybody \
	mc \
	mdi \
	meam \
	mesont \
	mgpt \
	misc \
	ml-hdnnp \
	ml-pace \
	ml-pod \
	ml-quip \
	ml-rann \
	ml-snap \
	ml-uf3 \
	mofff \
	molfile \
	netcdf \
	openmp \
	opt \
	orient \
	peri \
	plugin \
	plumed \
	poems \
	ptm \
	python \
	qeq \
	qmmm \
	qtb \
	reaction \
	reaxff \
	replica \
	rigid \
	scafacos \
	shock \
	smtbq \
	sph \
	spin \
	srd \
	tally \
	uef \
	voronoi \
	vtk \
	yaff \
	atc \
	dielectric \
	electrode \
	ml-iap \
	phonon

# NOTE: the last four packages must remain at the end since
# they depend on other packages to be installed first.

PACKBASIC = molecule kspace manybody rigid

PACKMOST = \
	kspace \
	molecule \
	amoeba \
	asphere \
	bocs \
	body \
	bpm \
	brownian \
	cg-dna \
	cg-spica \
	class2 \
	colloid \
	coreshell \
	diffraction \
	dipole \
	dpd-basic \
	dpd-meso \
	dpd-react \
	dpd-smooth \
	drude \
	eff \
	extra-command \
	extra-compute \
	extra-dump \
	extra-fix \
	extra-molecule \
	extra-pair \
	fep \
	granular \
	interlayer \
	manybody \
	mc \
	meam \
	mesont \
	misc \
	ml-snap \
	mofff \
	openmp \
	opt \
	orient \
	peri \
	plugin \
	qeq \
	reaction \
	reaxff \
	replica \
	rigid \
	shock \
	sph \
	spin \
	srd \
	uef \
	yaff \
	dielectric \
	ml-iap \
	phonon

# NOTE ^^^^^: the last three packages must remain at the end since
# they depend on other packages to be installed first.

PACKLIB = \
	compress \
	gpu \
	kim \
	kokkos \
	lepton \
	poems \
	python \
	voronoi \
	adios \
	atc \
	awpmd \
	colvars \
	h5md \
	ml-hdnnp \
	latboltz \
	lepton \
	mdi \
	molfile \
	netcdf \
	ml-pace \
	ml-pod \
	plumed \
	qmmm \
	ml-quip \
	scafacos \
	machdyn \
	vtk \
	electrode

PACKSYS = compress latboltz python

PACKINT = atc awpmd colvars electrode gpu kokkos lepton ml-pod poems

PACKEXT = \
	adios \
	h5md \
	kim \
	machdyn \
	ml-hdnnp \
	ml-pace \
	ml-quip \
	molfile \
	netcdf \
	plumed \
	qmmm \
	scafacos \
	voronoi \
	vtk \

# Helper GNU make function for conversion to upper case without using shell commands
uppercase_TABLE:=a,A b,B c,C d,D e,E f,F g,G h,H i,I j,J k,K l,L m,M n,N o,O p,P q,Q r,R s,S t,T u,U v,V w,W x,X y,Y z,Z
uppercase_internal=$(if $1,$$(subst $(firstword $1),$(call uppercase_internal,$(wordlist 2,$(words $1),$1),$2)),$2)
uppercase=$(eval uppercase_RESULT:=$(call uppercase_internal,$(uppercase_TABLE),$1))$(uppercase_RESULT)

PACKAGEUC = $(call uppercase,$(PACKAGE))
PACKAGESORTED = $(sort $(PACKAGEUC))

YESDIR = $(call uppercase,$(@:yes-%=%))
NODIR  = $(call uppercase,$(@:no-%=%))
LIBDIR = $(@:lib-%=%)

# List of all targets

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo 'make mpi-stubs           build dummy MPI library in STUBS'
	@echo 'make install-python      install LAMMPS wrapper in Python'
	@echo 'make tar                 create lmp_src.tar.gz for src dir and packages'
	@echo ''
	@echo 'make package                 list available packages and their dependencies'
	@echo 'make package-status (ps)     status of all packages'
	@echo 'make package-installed (pi)  list of installed packages'
	@echo 'make yes-package             install a single pkg in src dir'
	@echo 'make no-package              remove a single pkg from src dir'
	@echo 'make yes-all                 install all pkgs in src dir'
	@echo 'make no-all                  remove all pkgs from src dir'
	@echo 'make no-installed (ni)       quick uninstall of installed pkgs'
	@echo 'make yes-basic               install a few commonly used pkgs'
	@echo 'make no-basic                remove a few commonly used pkgs'
	@echo 'make yes-most                install most pkgs w/o libs'
	@echo 'make no-most                 remove most pkgs w/o libs'
	@echo 'make yes-lib       install all pkgs with libs (included or ext)'
	@echo 'make no-lib        remove all pkgs with libs (included or ext)'
	@echo 'make yes-ext                 install all pkgs with external libs'
	@echo 'make no-ext                  remove all pkgs with external libs'
	@echo ''
	@echo 'make package-update (pu) replace src files with updated package files'
	@echo 'make package-overwrite   replace package files with src files'
	@echo 'make package-diff (pd)   diff src files against package files'
	@echo ''
	@echo 'make lib-package         help for download/build/install a package library'
	@echo 'make lib-package args="..."    download/build/install a package library'
	@echo 'make purge               purge obsolete copies of source files'
	@echo ''
	@echo 'make machine             build LAMMPS for machine with static library'
	@echo 'make mode=static machine same as above'
	@echo 'make mode=shared machine build LAMMPS for machine with shared library'
	@echo 'make mode=print machine  print compiler/linker flags'
	@echo ''
	@echo 'machine is one of these from src/MAKE:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''
	@echo '... or one of these from src/MAKE/OPTIONS:'
	@echo ''
	@files="`ls MAKE/OPTIONS/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''
	@echo '... or one of these from src/MAKE/MACHINES:'
	@echo ''
	@files="`ls MAKE/MACHINES/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''
	@echo '... or one of these from src/MAKE/MINE:'
	@echo ''
	@files="`ls MAKE/MINE/Makefile.* 2>/dev/null`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''


lmpinstalledpkgs.h: $(SRC) $(INC)
	@echo '#ifndef LMP_INSTALLED_PKGS_H' >  ${TMPNAME}.lmpinstalled
	@echo '#define LMP_INSTALLED_PKGS_H' >> ${TMPNAME}.lmpinstalled
	@echo 'const char * LAMMPS_NS::LAMMPS::installed_packages[] = {' >> ${TMPNAME}.lmpinstalled
	@for p in $(PACKAGEUC); do info=$$($(SHELL) Package.sh $$p installed); \
             [ -n "$$info" ] && echo "\"$$info\"" | sed -e 's/".*package \(.*\)"/"\1",/' >> ${TMPNAME}.lmpinstalled || :; done
	@echo ' NULL };' >> ${TMPNAME}.lmpinstalled
	@echo '#endif' >> ${TMPNAME}.lmpinstalled
	@if [ -f lmpinstalledpkgs.h ]; \
          then test "`diff --brief ${TMPNAME}.lmpinstalled lmpinstalledpkgs.h`" != "" && \
	        mv ${TMPNAME}.lmpinstalled lmpinstalledpkgs.h || rm ${TMPNAME}.lmpinstalled ; \
        else mv ${TMPNAME}.lmpinstalled lmpinstalledpkgs.h ; fi

gitversion:
	@echo 'Gathering git version information'
	@echo '#ifndef LMP_GIT_VERSION_H' >  ${TMPNAME}.lmpgitversion
	@echo '#define LMP_GIT_VERSION_H' >> ${TMPNAME}.lmpgitversion
	@if (type git && test -e ../.git ) >> /dev/null 2>> /dev/null ; then \
	  git='true';					\
	  commit=$$(git rev-parse HEAD);		\
	  branch=$$(git rev-parse --abbrev-ref HEAD);	\
	  describe=$$(git describe --dirty=-modified);	\
	else \
	  git='false' ;					\
	  commit='(unknown)' ;				\
	  branch='(unknown)' ;				\
	  describe='(unknown)' ;			\
	fi ; \
	echo "bool LAMMPS_NS::LAMMPS::has_git_info() { return $${git}; }"        >> ${TMPNAME}.lmpgitversion ; \
	echo "const char *LAMMPS_NS::LAMMPS::git_commit() { return \"$${commit}\"; }" >> ${TMPNAME}.lmpgitversion ; \
	echo "const char *LAMMPS_NS::LAMMPS::git_branch() { return \"$${branch}\"; }" >> ${TMPNAME}.lmpgitversion ; \
	echo "const char *LAMMPS_NS::LAMMPS::git_descriptor() { return \"$${describe}\"; }" >> ${TMPNAME}.lmpgitversion
	@echo '#endif' >> ${TMPNAME}.lmpgitversion
	@if [ -f lmpgitversion.h ]; \
          then test "`diff --brief ${TMPNAME}.lmpgitversion lmpgitversion.h`" != "" && \
	        mv ${TMPNAME}.lmpgitversion lmpgitversion.h || rm ${TMPNAME}.lmpgitversion ; \
        else mv ${TMPNAME}.lmpgitversion lmpgitversion.h ; fi

# Build LAMMPS in one of 2 modes
# static = static compile in Obj_machine (default)
# shared = shared compile in Obj_shared_machine

.DEFAULT:
	@if [ $@ = "serial" ]; \
	  then cd STUBS; $(MAKE); cd ..; fi
	@test -f MAKE/Makefile.$@ -o -f MAKE/OPTIONS/Makefile.$@ -o \
	  -f MAKE/MACHINES/Makefile.$@ -o -f MAKE/MINE/Makefile.$@
	@if [ ! -d $(objdir) ]; then mkdir $(objdir); fi
	@echo 'Gathering installed package information (may take a little while)'
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh packages
	@$(MAKE) $(MFLAGS) lmpinstalledpkgs.h gitversion
	@echo 'Compiling LAMMPS for machine $@'
	@if [ -f MAKE/MACHINES/Makefile.$@ ]; \
	  then cp MAKE/MACHINES/Makefile.$@ $(objdir)/Makefile; fi
	@if [ -f MAKE/OPTIONS/Makefile.$@ ]; \
	  then cp MAKE/OPTIONS/Makefile.$@ $(objdir)/Makefile; fi
	@if [ -f MAKE/Makefile.$@ ]; \
	  then cp MAKE/Makefile.$@ $(objdir)/Makefile; fi
	@if [ -f MAKE/MINE/Makefile.$@ ]; \
	  then cp MAKE/MINE/Makefile.$@ $(objdir)/Makefile; fi
	@if [ ! -e Makefile.package ]; \
	  then cp Makefile.package.empty Makefile.package; fi
	@if [ ! -e Makefile.package.settings ]; \
	  then cp Makefile.package.settings.empty Makefile.package.settings; fi
	@cp Makefile.package Makefile.package.settings $(objdir)
	@cd $(objdir); rm -f .depend; \
	$(MAKE) $(MFLAGS) "SRC = $(SRC)" "INC = $(INC)" depend || :
	@rm -f $(ARLINK) $(SHLINK) $(EXE)
ifeq ($(mode),static)
	@cd $(objdir); \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJLIB)" "INC = $(INC)" "SHFLAGS =" \
	  "LMPLIB = $(ARLIB)" "ARLIB = $(ARLIB)" "SHLIB = $(SHLIB)" \
	  "LMPLINK = $(LMPLINK)" "DYN_LIB = $(DYN_LIB)" "EXE = ../$(EXE)" ../$(EXE)
	@ln -s $(ARLIB) $(ARLINK)
endif
ifeq ($(mode),shared)
	@cd $(objdir); \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJLIB)" "INC = $(INC)" \
	  "LMPLIB = $(SHLIB)" "ARLIB = $(ARLIB)" "SHLIB = $(SHLIB)" \
	  "LMPLINK = $(LMPLINK)" "DYN_LIB = $(DYN_LIB)" "EXE = ../$(EXE)" ../$(EXE)
	@ln -s $(SHLIB) $(SHLINK)
endif
# backward compatibility
ifeq ($(mode),exe)
	$(MAKE) $(MFLAGS) mode=static $@
endif
ifeq ($(mode),lib)
	$(MAKE) $(MFLAGS) mode=static $@
endif
ifeq ($(mode),shexe)
	$(MAKE) $(MFLAGS) mode=shared $@
endif
ifeq ($(mode),shlib)
	$(MAKE) $(MFLAGS) mode=shared $@
endif

ifeq ($(mode),print)
	@cd $(objdir); \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJLIB)" "INC = $(INC)" \
	  "EXE = ../$(ARLIB)" -f ../Makefile.print
endif

# Remove machine-specific object files

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -rf Obj_*
	rm -f style_*.h packages_*.h lmpgitversion.h lmpinstalledpkgs.h

clean-%:
	@if [ $@ = "clean-serial" ]; \
		then cd STUBS; $(MAKE) clean; cd ..; fi
	rm -rf Obj_$(@:clean-%=%) Obj_shared_$(@:clean-%=%)

# Make MPI STUBS library

mpi-stubs:
	@cd STUBS; $(MAKE) clean; $(MAKE)

# install LAMMPS shared lib and Python wrapper for Python usage
# include python package settings to automatically adapt name of
# the python interpreter. must purge build folder to not install
# unwanted outdated files.

sinclude ../lib/python/Makefile.lammps
install-python:
	@rm -rf ../python/build
	@$(PYTHON) ../python/install.py -p ../python/lammps -l ../src/liblammps.so -w $(PWD) -v $(PWD)/version.h

# Create a tarball of src dir and packages

tar:
	@cd STUBS; $(MAKE) clean
	@cd ..; tar cvzf src/$(ROOT)_src.tar.gz \
	  src/Make* src/Package.sh src/Depend.sh src/Install.sh src/Fetch.sh \
	  src/MAKE src/DEPEND src/*.cpp src/*.h src/STUBS \
	  $(patsubst %,src/%,$(PACKAGEUC)) \
          --exclude=*/.svn
	@cd STUBS; $(MAKE)
	@echo "Created $(ROOT)_src.tar.gz"

check: check-whitespace check-permissions check-homepage check-errordocs check-docs check-version

check-whitespace:
	$(PYTHON) ../tools/coding_standard/whitespace.py ..

fix-whitespace:
	$(PYTHON) ../tools/coding_standard/whitespace.py .. -f

check-permissions:
	$(PYTHON) ../tools/coding_standard/permissions.py ..

fix-permissions:
	$(PYTHON) ../tools/coding_standard/permissions.py .. -f

check-homepage:
	$(PYTHON) ../tools/coding_standard/homepage.py ..

fix-homepage:
	$(PYTHON) ../tools/coding_standard/homepage.py .. -f

check-errordocs:
	$(PYTHON) ../tools/coding_standard/errordocs.py ..

fix-errordocs:
	$(PYTHON) ../tools/coding_standard/errordocs.py .. -f

check-docs:
	$(MAKE) $(MFLAGS) -C ../doc anchor_check style_check package_check role_check

check-version:
	$(PYTHON) ../tools/coding_standard/versiontags.py .. || echo


format-src:
	clang-format -i --verbose --style=file *.cpp *.h */*.cpp */*.h

format-tests:
	clang-format -i --verbose --style=file ../unittest/*/*.cpp ../unittest/*/*.h


# Package management

package:
	@echo 'Available packages:' $(PACKAGE)
	@echo ''
	@echo 'Packages that need system libraries:' $(PACKSYS)
	@echo ''
	@echo 'Packages that need provided libraries:' $(PACKINT)
	@echo ''
	@echo 'Packages that need external libraries:' $(PACKEXT)
	@echo ''
	@echo 'make package                 list available packages'
	@echo 'make package                 list available packages'
	@echo 'make package-status (ps)     status of all packages'
	@echo 'make package-installed (pi)  list of installed packages'
	@echo 'make yes-package             install a single pgk in src dir'
	@echo 'make no-package              remove a single pkg from src dir'
	@echo 'make yes-basic               install a few commonly used pgks in src'
	@echo 'make no-basic                remove a few commonly used pkgs from src dir'
	@echo 'make yes-most                install most pgks w/o libs in src'
	@echo 'make no-most                 remove most pkgs w/o libs from src dir'
	@echo 'make yes-all                 install all pgks in src dir'
	@echo 'make no-all                  remove all pkgs from src dir'
	@echo 'make yes-lib                 install all pkgs with libs (included or ext)'
	@echo 'make no-lib                  remove all pkgs with libs (included or ext)'
	@echo 'make yes-ext                 install all pkgs with external libs'
	@echo 'make no-ext                  remove all pkgs with external libs'
	@echo ''
	@echo 'make package-update (pu)  replace src files with package files'
	@echo 'make package-overwrite    replace package files with src files'
	@echo 'make package-diff (pd)    diff src files against package file'
	@echo ''
	@echo 'make lib-package      build and/or download a package library'

yes-all:
	@for p in $(PACKAGE); do $(MAKE) yes-$$p; done

no-all:
	@for p in $(PACKAGE); do $(MAKE) no-$$p; done

no-installed ni:
	@for p in $(PACKAGESORTED); do $(SHELL) Package.sh $$p info && $(MAKE) no-$$p || : ; done

yes-standard yes-std:
	@echo 'There are no more "standard" or "user" packages in LAMMPS'

no-standard no-std:
	@echo 'There are no more "standard" or "user" packages in LAMMPS'

yes-user:
	@echo 'There are no more "standard" or "user" packages in LAMMPS'

no-user:
	@echo 'There are no more "standard" or "user" packages in LAMMPS'

yes-basic:
	@for p in $(PACKBASIC); do $(MAKE) yes-$$p; done

no-basic:
	@for p in $(PACKBASIC); do $(MAKE) no-$$p; done

yes-most:
	@for p in $(PACKMOST); do $(MAKE) yes-$$p; done

no-most:
	@for p in $(PACKMOST); do $(MAKE) no-$$p; done

yes-lib:
	@for p in $(PACKLIB); do $(MAKE) yes-$$p; done

no-lib:
	@for p in $(PACKLIB); do $(MAKE) no-$$p; done

yes-ext:
	@for p in $(PACKEXT); do $(MAKE) yes-$$p; done

no-ext:
	@for p in $(PACKEXT); do $(MAKE) no-$$p; done

yes-%:
	@if [ ! -e Makefile.package ]; \
	  then cp Makefile.package.empty Makefile.package; fi
	@if [ ! -e Makefile.package.settings ]; \
	  then cp Makefile.package.settings.empty Makefile.package.settings; fi
	@if [ ! -e $(YESDIR) ]; then \
	  echo "Package $(YESDIR) does not exist"; exit 1; \
	elif [ -e $(YESDIR)/Install.sh ]; then \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); $(SHELL) Install.sh 1; cd ..; \
		$(SHELL) Depend.sh $(YESDIR) 1; \
	else \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); $(SHELL) ../Install.sh 1; cd ..; \
		$(SHELL) Depend.sh $(YESDIR) 1; \
	fi;
	@$(SHELL) Fetch.sh $(YESDIR)

no-%:
	@if [ ! -e $(NODIR) ]; then \
	  echo "Package $(NODIR) does not exist"; exit 1; \
	elif [ -e $(NODIR)/Install.sh ]; then \
	  echo "Uninstalling package $(@:no-%=%)"; \
	  cd $(NODIR); $(SHELL) Install.sh 0; cd ..; \
		$(SHELL) Depend.sh $(NODIR) 0; \
	else \
	  echo "Uninstalling package $(@:no-%=%)"; \
	  cd $(NODIR); $(SHELL) ../Install.sh 0; cd ..; \
		$(SHELL) Depend.sh $(NODIR) 0; \
        fi;

# download/build/install a package library
# update the timestamp on main.cpp to trigger a relink with "make machine"

lib-%:
	@if [ -e ../lib/$(LIBDIR)/Install.py ]; then \
	  echo "Installing lib $(@:lib-%=%)"; \
	  ( cd ../lib/$(LIBDIR); $(PYTHON) Install.py $(args) ); \
	else \
	  echo "Install script for lib $(@:lib-%=%) does not exist"; \
	fi; touch main.cpp

lib-colvars: lib-lepton

# status = list src files that differ from package files
# installed = list of installed packages
# update = replace src files with newer package files
# overwrite = overwrite package files with newer src files
# diff = show differences between src and package files
# purge = delete obsolete and auto-generated package files

package-status ps:
	@for p in $(PACKAGESORTED); do $(SHELL) Package.sh $$p status; done

package-installed pi:
	@for p in $(PACKAGESORTED); do $(SHELL) Package.sh $$p installed; done

package-update pu: purge
	@echo 'Updating installed packages:'
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p update; done

package-overwrite: purge
	@echo 'Overwriting installed packages:'
	@for p in $(PACKAGESORTED); do $(SHELL) Package.sh $$p overwrite; done

package-diff pd:
	@for p in $(PACKAGESORTED); do $(SHELL) Package.sh $$p diff; done

purge: Purge.list
	@echo 'Purging obsolete and auto-generated source files'
	@for f in `grep -v '#' Purge.list` ;		\
	    do test -f $$f && rm $$f && echo $$f || : ;		\
	done
