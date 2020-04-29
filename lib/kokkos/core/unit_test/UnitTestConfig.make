KOKKOS_PATH = ../..

# See $(KOKKOS_PATH)/Makefile.kokkos and $(KOKKOS_PATH)/generate_makefile.bash
KOKKOS_ARCH_OPTIONS="None AMDAVX ARMv80 ARMv81 ARMv8-ThunderX \
	 BGQ Power7 Power8 Power9 \
	 WSM SNB HSW BDW SKX KNC KNL \
     Kepler Kepler30 Kepler32 Kepler35 Kepler37 \
     Maxwell Maxwell50 Maxwell52 Maxwell53 Pascal60 Pascal61"
#KOKKOS_ARCH_OPTIONS="AMDAVX"

KOKKOS_DEVICE_OPTIONS="Cuda ROCm OpenMP Pthread Serial"
#KOKKOS_DEVICE_OPTIONS="Cuda"

# Configure paths to enable environment query in Makefile.kokkos to work
ROCM_HCC_PATH="config"
CXX="./config/cxx"
ipath=env CXX=$(CXX) env PATH=./config:$$PATH env ROCM_HCC_PATH=$(ROCM_HCC_PATH)

# Defined in core/src/Makefile -- this should be consistent
KOKKOS_MAKEFILE=Makefile.kokkos
KOKKOS_CMAKEFILE=kokkos_generated_settings.cmake

# Defined in Makefile.kokkos -- this should be consistent
KOKKOS_INTERNAL_CONFIG_TMP=KokkosCore_config.tmp
KOKKOS_CONFIG_HEADER=KokkosCore_config.h

d='\#'

# diff => 0 is no difference.  if => 0 is false
testmake=if test "`testmake.sh $1 $2 $3`" = 'Passed'; then echo OK $d $1; else echo not OK $d $1; fi
testconf=if test "`diffconfig.sh $1`" = 'Passed'; then echo OK $d $1; else echo not OK $d $1; fi

# testing tmp and cmakefile files is unnecessary here
test:
	@for karch in "$(KOKKOS_ARCH_OPTIONS)"; do \
	  for device in "$(KOKKOS_DEVICE_OPTIONS)"; do \
	     $(ipath) KOKKOS_DEVICES=$$device KOKKOS_ARCH=$$karch make -e -f ../src/Makefile build-makefile-cmake-kokkos; \
		 rm -f $(KOKKOS_INTERNAL_CONFIG_TMP) $(KOKKOS_CMAKEFILE); \
		 prfx="$$karch"_"$$device"_; \
		 newmake="$$prfx"$(KOKKOS_MAKEFILE);  \
		 newconf="$$prfx"$(KOKKOS_CONFIG_HEADER); \
		 mv $(KOKKOS_MAKEFILE)      config/tmpstore/$$newmake; \
		 mv $(KOKKOS_CONFIG_HEADER) config/tmpstore/$$newconf; \
		 $(call testmake,$$newmake,$$karch,$$device); \
		 $(call testconf,$$newconf); \
	  done; \
	done

test-cmake:
	@cd config/cmaketest; \
     cmake .             ; \
     make test
