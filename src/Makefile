# LAMMPS multiple-machine Makefile

SHELL = /bin/bash
#.IGNORE:

# Definitions

ROOT =	lmp
EXE =	$(ROOT)_$@
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# Package variables

PACKAGE = asphere body class2 colloid dipole fld gpu granular kim \
	  kokkos kspace manybody mc meam misc molecule mpiio opt peri poems \
	  qeq reax replica rigid shock snap srd voronoi xtc

PACKUSER = user-atc user-awpmd user-cg-cmm user-colvars \
	   user-cuda user-eff user-fep user-intel user-lb user-misc \
	   user-molfile user-omp user-phonon user-pimd \
	   user-qmmm user-reaxc user-sph

PACKLIB = gpu kim meam poems reax voronoi \
	  user-atc user-awpmd user-colvars user-pimd \
	  user-qmmm user-cuda user-molfile

PACKALL = $(PACKAGE) $(PACKUSER)

PACKAGEUC = $(shell echo $(PACKAGE) | tr a-z A-Z)
PACKUSERUC = $(shell echo $(PACKUSER) | tr a-z A-Z)

YESDIR = $(shell echo $(@:yes-%=%) | tr a-z A-Z)
NODIR  = $(shell echo $(@:no-%=%) | tr a-z A-Z)

# List of all targets

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo 'make makelib             create Makefile.lib for static library build'
	@echo 'make makeshlib           create Makefile.shlib for shared library build'
	@echo 'make makelist            create Makefile.list used by old makes'
	@echo 'make -f Makefile.lib machine      build LAMMPS as static library for machine'
	@echo 'make -f Makefile.shlib machine    build LAMMPS as shared library for machine'
	@echo 'make -f Makefile.list machine     build LAMMPS from explicit list of files'
	@echo 'make stubs               build dummy MPI library in STUBS'
	@echo 'make install-python      install LAMMPS wrapper in Python'
	@echo 'make tar                 create lmp_src.tar.gz of src dir and packages'
	@echo ''
	@echo 'make package                 list available packages'
	@echo 'make package-status (ps)     status of all packages'
	@echo 'make yes-package             install a single pgk in src dir'
	@echo 'make no-package              remove a single pkg from src dir'
	@echo 'make yes-all                 install all pgks in src dir'
	@echo 'make no-all                  remove all pkgs from src dir'
	@echo 'make yes-standard (yes-std)  install all standard pkgs'
	@echo 'make no-standard (no-std)    remove all standard pkgs'
	@echo 'make yes-user                install all user pkgs'
	@echo 'make no-user                 remove all user pkgs'
	@echo 'make no-lib                  remove all pkgs with external libs'
	@echo ''
	@echo 'make package-update (pu) replace src files with updated package files'
	@echo 'make package-overwrite   replace package files with src files'
	@echo 'make package-diff (pd)   diff src files against package files'
	@echo 'make package-purge       purge obsolete copies of package sources'
	@echo ''
	@echo 'make machine             build LAMMPS where machine is ...'
	@echo ''
	@echo 'one of these from src/MAKE:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''
	@echo 'or one of these from src/MAKE/OPTIONS:'
	@echo ''
	@files="`ls MAKE/OPTIONS/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''
	@echo 'or one of these from src/MAKE/MACHINES:'
	@echo ''
	@files="`ls MAKE/MACHINES/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''
	@echo 'or one of these from src/MAKE/MINE:'
	@echo ''
	@files="`ls MAKE/MINE/Makefile.* 2>/dev/null`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''

# Build LAMMPS

.DEFAULT:
	@if [ $@ = "serial" -a ! -f STUBS/libmpi_stubs.a ]; \
	  then $(MAKE) stubs; fi
	@test -f MAKE/Makefile.$@ -o -f MAKE/OPTIONS/Makefile.$@ -o \
	  -f MAKE/MACHINES/Makefile.$@ -o -f MAKE/MINE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@$(SHELL) Make.sh style
	@if [ -f MAKE/MACHINES/Makefile.$@ ]; \
	  then cp MAKE/MACHINES/Makefile.$@ Obj_$@/Makefile; fi
	@if [ -f MAKE/OPTIONS/Makefile.$@ ]; \
	  then cp MAKE/OPTIONS/Makefile.$@ Obj_$@/Makefile; fi
	@if [ -f MAKE/Makefile.$@ ]; \
	  then cp MAKE/Makefile.$@ Obj_$@/Makefile; fi
	@if [ -f MAKE/MINE/Makefile.$@ ]; \
	  then cp MAKE/MINE/Makefile.$@ Obj_$@/Makefile; fi
	@if [ ! -e Makefile.package ]; \
	  then cp Makefile.package.empty Makefile.package; fi
	@if [ ! -e Makefile.package.settings ]; \
	  then cp Makefile.package.settings.empty Makefile.package.settings; fi
	@cp Makefile.package Makefile.package.settings Obj_$@
	@cd Obj_$@; \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" "SHFLAGS =" \
	  "EXE = ../$(EXE)" ../$(EXE)

# Remove machine-specific object files

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -rf Obj_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)

# Create Makefile.lib, Makefile.shlib, and Makefile.list

makelib:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.lib

makeshlib:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.shlib

makelist:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.list

# Make MPI STUBS library

stubs:
	@cd STUBS; $(MAKE) clean; $(MAKE)

# install LAMMPS shared lib and Python wrapper for Python usage

install-python:
	@python ../python/install.py

# Create a tarball of src dir and packages

tar:
	@cd STUBS; $(MAKE) clean
	@cd ..; tar cvzf src/$(ROOT)_src.tar.gz \
	  src/Make* src/Package.sh src/MAKE src/*.cpp src/*.h src/STUBS \
	  $(patsubst %,src/%,$(PACKAGEUC)) $(patsubst %,src/%,$(PACKUSERUC)) \
          --exclude=*/.svn
	@cd STUBS; $(MAKE)
	@echo "Created $(ROOT)_src.tar.gz"

# Package management

package:
	@echo 'Standard packages:' $(PACKAGE)
	@echo ''
	@echo 'User-contributed packages:' $(PACKUSER)
	@echo ''
	@echo 'make package                 list available packages'
	@echo 'make package-status (ps)     status of all packages'
	@echo 'make yes-package             install a single pgk in src dir'
	@echo 'make no-package              remove a single pkg from src dir'
	@echo 'make yes-all                 install all pgks in src dir'
	@echo 'make no-all                  remove all pkgs from src dir'
	@echo 'make yes-standard (yes-std)  install all standard pkgs'
	@echo 'make no-standard (no-srd)    remove all standard pkgs'
	@echo ''
	@echo 'make yes-user                install all user pkgs'
	@echo 'make no-user                 remove all user pkgs'
	@echo 'make no-lib                  remove all pkgs with external libs'
	@echo 'make package-update (pu)  replace src files with package files'
	@echo 'make package-overwrite    replace package files with src files'
	@echo 'make package-diff (pd)    diff src files against package file'

yes-all:
	@for p in $(PACKALL); do $(MAKE) yes-$$p; done

no-all:
	@for p in $(PACKALL); do $(MAKE) no-$$p; done

yes-standard yes-std:
	@for p in $(PACKAGE); do $(MAKE) yes-$$p; done

no-standard no-std:
	@for p in $(PACKAGE); do $(MAKE) no-$$p; done

yes-user:
	@for p in $(PACKUSER); do $(MAKE) yes-$$p; done

no-user:
	@for p in $(PACKUSER); do $(MAKE) no-$$p; done

no-lib:
	@for p in $(PACKLIB); do $(MAKE) no-$$p; done

yes-%:
	@if [ ! -e Makefile.package ]; \
	  then cp Makefile.package.empty Makefile.package; fi
	@if [ ! -e Makefile.package.settings ]; \
	  then cp Makefile.package.settings.empty Makefile.package.settings; fi
	@if [ ! -e $(YESDIR) ]; then \
	  echo "Package $(@:yes-%=%) does not exist"; \
	elif [ -e $(YESDIR)/Install.sh ]; then \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); $(SHELL) Install.sh 1; cd ..; \
		$(SHELL) Depend.sh $(YESDIR) 1; \
	else \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); $(SHELL) ../Install.sh 1; cd ..; \
		$(SHELL) Depend.sh $(YESDIR) 1; \
	fi;

no-%:
	@if [ ! -e $(NODIR) ]; then \
	  echo "Package $(@:no-%=%) does not exist"; \
	elif [ -e $(NODIR)/Install.sh ]; then \
	  echo "Uninstalling package $(@:no-%=%)"; \
	  cd $(NODIR); $(SHELL) Install.sh 0; cd ..; \
		$(SHELL) Depend.sh $(NODIR) 0; \
	else \
	  echo "Uninstalling package $(@:no-%=%)"; \
	  cd $(NODIR); $(SHELL) ../Install.sh 0; cd ..; \
		$(SHELL) Depend.sh $(NODIR) 0; \
        fi;

# status = list src files that differ from package files
# update = replace src files with newer package files
# overwrite = overwrite package files with newer src files
# diff = show differences between src and package files
# purge = delete obsolete and auto-generated package files

package-status ps:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p status; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p status; done

package-update pu:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p update; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p update; done

package-overwrite:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p overwrite; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p overwrite; done

package-diff pd:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p diff; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p diff; done

package-purge: Purge.list
	@echo 'Purging obsolete and auto-generated source files'
	@for f in `grep -v '#' Purge.list` ;		\
	    do test -f $$f && rm $$f && echo $$f || : ;		\
	done
