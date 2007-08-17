# LAMMPS multiple-machine Makefile

SHELL = /bin/sh
#.IGNORE:

# Definitions

ROOT =	lmp
EXE =	$(ROOT)_$@
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# Package variables

PACKAGE = asphere class2 colloid dipole dpd granular \
	  kspace manybody meam molecule opt poems xtc
PACKUC = $(shell perl -e 'printf("%s", uc("$(PACKAGE)"));')
YESDIR = $(shell perl -e 'printf("%s", uc("$(@:yes-%=%)"));')
NODIR  = $(shell perl -e 'printf("%s", uc("$(@:no-%=%)"));')

# List of all targets

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo 'make tar                 lmp_src.tar.gz of src dir and packages'
	@echo 'make makelib             update Makefile.lib for library build'
	@echo 'make makelist            update Makefile.list used by old makes'
	@echo ''
	@echo 'make package             list available packages'
	@echo 'make package-status      status of all packages'
	@echo 'make yes-package         install a package in src dir'
	@echo 'make no-package          remove package files from src dir'
	@echo 'make yes-all             install all packages in src dir'
	@echo 'make no-all              remove all package files from src dir'
	@echo 'make package-update      replace src files with package files'
	@echo 'make package-overwrite   replace package files with src files'
	@echo ''
	@echo 'make machine             build LAMMPS where machine is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''

# Build the code

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp -p $(SRC) $(INC) Obj_$@
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE)  "OBJ = $(OBJ)" "INC = $(INC)" "EXE = ../$(EXE)" ../$(EXE)
	@if [ -d Obj_$@ ]; then cd Obj_$@; rm $(SRC) $(INC) Makefile*; fi

# Remove machine-specific object files

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -r Obj_*

clean-%:
	rm -r Obj_$(@:clean-%=%)

# Create a tarball of src dir and packages

tar:
	@cd STUBS; make clean
	@cd ..; tar cvzf src/$(ROOT)_src.tar.gz \
	  src/Make* src/Package.csh src/MAKE src/*.cpp src/*.h src/STUBS \
	  $(patsubst %,src/%,$(PACKUC)) --exclude=*/.svn
	@cd STUBS; make
	@echo "Created $(ROOT)_src.tar.gz"

# Update Makefile.lib and Makefile.list

makelib:
	@csh Make.csh Makefile.lib

makelist:
	@csh Make.csh Makefile.list

# Package management
# status =    list differences between src and package files
# update =    replace src files with newer package files
# overwrite = overwrite package files with newer src files

package:
	@echo 'Available packages:'
	@echo $(PACKAGE)
	@echo ''
	@echo 'make package-status      status of all packages'
	@echo 'make yes-package         install a package in src dir'
	@echo 'make no-package          remove package files from src dir'
	@echo 'make yes-all             install all packages in src dir'
	@echo 'make no-all              remove all package files from src dir'
	@echo 'make package-update      replace src files with package files'
	@echo 'make package-overwrite   replace package files with src files'

yes-all:
	@for p in $(PACKAGE); do $(MAKE) yes-$$p; done

no-all:
	@for p in $(PACKAGE); do $(MAKE) no-$$p; done

yes-%:
	@if [ ! -e $(YESDIR) ]; then \
	  echo "Package $(@:yes-%=%) does not exist"; \
	else \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); csh -f Install.csh 1; \
	fi;

no-%:
	@if [ ! -e $(NODIR) ]; then \
	  echo "Package $(@:no-%=%) does not exist"; \
	else \
	  echo "Uninstalling package $(@:no-%=%), ignore errors"; \
	  cd $(NODIR); csh -f Install.csh 0; cd ..; $(MAKE) clean-all; \
        fi;

package-status:
	@for p in $(PACKUC); do csh -f Package.csh $$p status; done

package-update:
	@for p in $(PACKUC); do csh -f Package.csh $$p update; done

package-overwrite:
	@for p in $(PACKUC); do csh -f Package.csh $$p overwrite; done
