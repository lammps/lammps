ifndef KOKKOS_PATH
  MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
  KOKKOS_PATH = $(subst Makefile,,$(MAKEFILE_PATH))../..
endif

PREFIX ?= /usr/local/lib/kokkos

default: messages build-lib
	echo "End Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
  CXX = $(KOKKOS_PATH)/bin/nvcc_wrapper
else
  CXX = g++
endif

CXXFLAGS = -O3
LINK ?= $(CXX)
LDFLAGS ?=

include $(KOKKOS_PATH)/Makefile.kokkos

PWD = $(shell pwd)

KOKKOS_HEADERS_INCLUDE = $(wildcard $(KOKKOS_PATH)/core/src/*.hpp)
KOKKOS_HEADERS_INCLUDE_IMPL = $(wildcard $(KOKKOS_PATH)/core/src/impl/*.hpp)
KOKKOS_HEADERS_INCLUDE += $(wildcard $(KOKKOS_PATH)/containers/src/*.hpp)
KOKKOS_HEADERS_INCLUDE_IMPL += $(wildcard $(KOKKOS_PATH)/containers/src/impl/*.hpp)
KOKKOS_HEADERS_INCLUDE += $(wildcard $(KOKKOS_PATH)/algorithms/src/*.hpp)

CONDITIONAL_COPIES =

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
  KOKKOS_HEADERS_CUDA += $(wildcard $(KOKKOS_PATH)/core/src/Cuda/*.hpp)
  CONDITIONAL_COPIES += copy-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)
  KOKKOS_HEADERS_THREADS += $(wildcard $(KOKKOS_PATH)/core/src/Threads/*.hpp)
  CONDITIONAL_COPIES += copy-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_QTHREADS), 1)
  KOKKOS_HEADERS_QTHREADS += $(wildcard $(KOKKOS_PATH)/core/src/Qthreads/*.hpp)
  CONDITIONAL_COPIES += copy-qthreads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
  KOKKOS_HEADERS_OPENMP += $(wildcard $(KOKKOS_PATH)/core/src/OpenMP/*.hpp)
  CONDITIONAL_COPIES += copy-openmp
endif

ifeq ($(KOKKOS_OS),CYGWIN)
  COPY_FLAG = -u
endif
ifeq ($(KOKKOS_OS),Linux)
  COPY_FLAG = -u
endif
ifeq ($(KOKKOS_OS),Darwin)
  COPY_FLAG =
endif

ifeq ($(KOKKOS_DEBUG),"no")
  KOKKOS_DEBUG_CMAKE = OFF
else
  KOKKOS_DEBUG_CMAKE = ON
endif

messages:
	echo "Start Build"

build-makefile-kokkos:
	rm -f Makefile.kokkos
	echo "#Global Settings used to generate this library" >> Makefile.kokkos
	echo "KOKKOS_PATH = $(PREFIX)" >> Makefile.kokkos
	echo "KOKKOS_DEVICES = $(KOKKOS_DEVICES)" >> Makefile.kokkos
	echo "KOKKOS_ARCH = $(KOKKOS_ARCH)" >> Makefile.kokkos
	echo "KOKKOS_DEBUG = $(KOKKOS_DEBUG)" >> Makefile.kokkos
	echo "KOKKOS_USE_TPLS = $(KOKKOS_USE_TPLS)" >> Makefile.kokkos
	echo "KOKKOS_CXX_STANDARD = $(KOKKOS_CXX_STANDARD)" >> Makefile.kokkos
	echo "KOKKOS_OPTIONS = $(KOKKOS_OPTIONS)" >> Makefile.kokkos
	echo "KOKKOS_CUDA_OPTIONS = $(KOKKOS_CUDA_OPTIONS)" >> Makefile.kokkos
	echo "CXX ?= $(CXX)" >> Makefile.kokkos
	echo "NVCC_WRAPPER ?= $(PREFIX)/bin/nvcc_wrapper" >> Makefile.kokkos
	echo "" >> Makefile.kokkos
	echo "#Source and Header files of Kokkos relative to KOKKOS_PATH" >> Makefile.kokkos
	echo "KOKKOS_HEADERS = $(KOKKOS_HEADERS)" >> Makefile.kokkos
	echo "KOKKOS_SRC = $(KOKKOS_SRC)" >> Makefile.kokkos
	echo "" >> Makefile.kokkos
	echo "#Variables used in application Makefiles" >> Makefile.kokkos
	echo "KOKKOS_OS = $(KOKKOS_OS)" >> Makefile.kokkos
	echo "KOKKOS_CPP_DEPENDS = $(KOKKOS_CPP_DEPENDS)" >> Makefile.kokkos
	echo "KOKKOS_CXXFLAGS = $(KOKKOS_CXXFLAGS)" >> Makefile.kokkos
	echo "KOKKOS_CPPFLAGS = $(KOKKOS_CPPFLAGS)" >> Makefile.kokkos
	echo "KOKKOS_LINK_DEPENDS  = $(KOKKOS_LINK_DEPENDS)" >> Makefile.kokkos
	echo "KOKKOS_LIBS = $(KOKKOS_LIBS)" >> Makefile.kokkos
	echo "KOKKOS_LDFLAGS = $(KOKKOS_LDFLAGS)" >> Makefile.kokkos
	echo "" >> Makefile.kokkos
	echo "#Internal settings which need to propagated for Kokkos examples" >> Makefile.kokkos
	echo "KOKKOS_INTERNAL_USE_CUDA = ${KOKKOS_INTERNAL_USE_CUDA}" >> Makefile.kokkos
	echo "KOKKOS_INTERNAL_USE_QTHREADS = ${KOKKOS_INTERNAL_USE_QTHREADS}" >> Makefile.kokkos
	echo "KOKKOS_INTERNAL_USE_OPENMP = ${KOKKOS_INTERNAL_USE_OPENMP}" >> Makefile.kokkos
	echo "KOKKOS_INTERNAL_USE_PTHREADS = ${KOKKOS_INTERNAL_USE_PTHREADS}" >> Makefile.kokkos
	echo "" >> Makefile.kokkos
	echo "#Fake kokkos-clean target" >> Makefile.kokkos
	echo "kokkos-clean:" >> Makefile.kokkos
	echo "" >> Makefile.kokkos
	sed \
		-e 's|$(KOKKOS_PATH)/core/src|$(PREFIX)/include|g' \
		-e 's|$(KOKKOS_PATH)/containers/src|$(PREFIX)/include|g' \
		-e 's|$(KOKKOS_PATH)/algorithms/src|$(PREFIX)/include|g' \
		-e 's|-L$(PWD)|-L$(PREFIX)/lib|g' \
		-e 's|= libkokkos.a|= $(PREFIX)/lib/libkokkos.a|g' \
		-e 's|= KokkosCore_config.h|= $(PREFIX)/include/KokkosCore_config.h|g' Makefile.kokkos \
		> Makefile.kokkos.tmp
	mv -f Makefile.kokkos.tmp Makefile.kokkos

build-cmake-kokkos:
	rm -f kokkos.cmake
	echo "#Global Settings used to generate this library" >> kokkos.cmake
	echo "set(KOKKOS_PATH $(PREFIX) CACHE PATH \"Kokkos installation path\")" >> kokkos.cmake
	echo "set(KOKKOS_DEVICES $(KOKKOS_DEVICES) CACHE STRING \"Kokkos devices list\")" >> kokkos.cmake
	echo "set(KOKKOS_ARCH $(KOKKOS_ARCH) CACHE STRING \"Kokkos architecture flags\")" >> kokkos.cmake
	echo "set(KOKKOS_DEBUG $(KOKKOS_DEBUG_CMAKE) CACHE BOOL \"Kokkos debug enabled ?)\")" >> kokkos.cmake
	echo "set(KOKKOS_USE_TPLS $(KOKKOS_USE_TPLS) CACHE STRING \"Kokkos templates list\")" >> kokkos.cmake
	echo "set(KOKKOS_CXX_STANDARD $(KOKKOS_CXX_STANDARD) CACHE STRING \"Kokkos C++ standard\")" >> kokkos.cmake
	echo "set(KOKKOS_OPTIONS $(KOKKOS_OPTIONS) CACHE STRING \"Kokkos options\")" >> kokkos.cmake
	echo "set(KOKKOS_CUDA_OPTIONS $(KOKKOS_CUDA_OPTIONS) CACHE STRING \"Kokkos Cuda options\")" >> kokkos.cmake
	echo "if(NOT $ENV{CXX})" >> kokkos.cmake
	echo '  message(WARNING "You are currently using compiler $${CMAKE_CXX_COMPILER} while Kokkos was built with $(CXX) ; make sure this is the behavior you intended to be.")' >> kokkos.cmake
	echo "endif()" >> kokkos.cmake
	echo "if(NOT DEFINED ENV{NVCC_WRAPPER})" >> kokkos.cmake
	echo "  set(NVCC_WRAPPER \"$(NVCC_WRAPPER)\" CACHE FILEPATH \"Path to command nvcc_wrapper\")" >> kokkos.cmake
	echo "else()" >> kokkos.cmake
	echo '  set(NVCC_WRAPPER $$ENV{NVCC_WRAPPER} CACHE FILEPATH "Path to command nvcc_wrapper")' >> kokkos.cmake
	echo "endif()" >> kokkos.cmake
	echo "" >> kokkos.cmake
	echo "#Source and Header files of Kokkos relative to KOKKOS_PATH" >> kokkos.cmake
	echo "set(KOKKOS_HEADERS \"$(KOKKOS_HEADERS)\" CACHE STRING \"Kokkos headers list\")" >> kokkos.cmake
	echo "set(KOKKOS_SRC \"$(KOKKOS_SRC)\" CACHE STRING \"Kokkos source list\")" >> kokkos.cmake
	echo "" >> kokkos.cmake
	echo "#Variables used in application Makefiles" >> kokkos.cmake
	echo "set(KOKKOS_CPP_DEPENDS \"$(KOKKOS_CPP_DEPENDS)\" CACHE STRING \"\")" >> kokkos.cmake
	echo "set(KOKKOS_CXXFLAGS \"$(KOKKOS_CXXFLAGS)\" CACHE STRING \"\")" >> kokkos.cmake
	echo "set(KOKKOS_CPPFLAGS \"$(KOKKOS_CPPFLAGS)\" CACHE STRING \"\")" >> kokkos.cmake
	echo "set(KOKKOS_LINK_DEPENDS \"$(KOKKOS_LINK_DEPENDS)\" CACHE STRING \"\")" >> kokkos.cmake
	echo "set(KOKKOS_LIBS \"$(KOKKOS_LIBS)\" CACHE STRING \"\")" >> kokkos.cmake
	echo "set(KOKKOS_LDFLAGS \"$(KOKKOS_LDFLAGS)\" CACHE STRING \"\")" >> kokkos.cmake
	echo "" >> kokkos.cmake
	echo "#Internal settings which need to propagated for Kokkos examples" >> kokkos.cmake
	echo "set(KOKKOS_INTERNAL_USE_CUDA \"${KOKKOS_INTERNAL_USE_CUDA}\" CACHE STRING \"\")" >> kokkos.cmake
	echo "set(KOKKOS_INTERNAL_USE_OPENMP \"${KOKKOS_INTERNAL_USE_OPENMP}\" CACHE STRING \"\")" >> kokkos.cmake
	echo "set(KOKKOS_INTERNAL_USE_PTHREADS \"${KOKKOS_INTERNAL_USE_PTHREADS}\" CACHE STRING \"\")" >> kokkos.cmake
	echo "mark_as_advanced(KOKKOS_HEADERS KOKKOS_SRC KOKKOS_INTERNAL_USE_CUDA KOKKOS_INTERNAL_USE_OPENMP KOKKOS_INTERNAL_USE_PTHREADS)" >> kokkos.cmake
	echo "" >> kokkos.cmake
	sed \
		-e 's|$(KOKKOS_PATH)/core/src|$(PREFIX)/include|g' \
	 	-e 's|$(KOKKOS_PATH)/containers/src|$(PREFIX)/include|g' \
	 	-e 's|$(KOKKOS_PATH)/algorithms/src|$(PREFIX)/include|g' \
	 	-e 's|-L$(PWD)|-L$(PREFIX)/lib|g' \
	 	-e 's|= libkokkos.a|= $(PREFIX)/lib/libkokkos.a|g' \
	 	-e 's|= KokkosCore_config.h|= $(PREFIX)/include/KokkosCore_config.h|g' kokkos.cmake \
	 	> kokkos.cmake.tmp
	mv -f kokkos.cmake.tmp kokkos.cmake

build-lib: build-makefile-kokkos build-cmake-kokkos $(KOKKOS_LINK_DEPENDS)

mkdir:
	mkdir -p $(PREFIX)
	mkdir -p $(PREFIX)/bin
	mkdir -p $(PREFIX)/include
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include/impl

copy-cuda: mkdir
	mkdir -p $(PREFIX)/include/Cuda
	cp $(COPY_FLAG) $(KOKKOS_HEADERS_CUDA) $(PREFIX)/include/Cuda

copy-threads: mkdir
	mkdir -p $(PREFIX)/include/Threads
	cp $(COPY_FLAG) $(KOKKOS_HEADERS_THREADS) $(PREFIX)/include/Threads

copy-qthreads: mkdir
	mkdir -p $(PREFIX)/include/Qthreads
	cp $(COPY_FLAG) $(KOKKOS_HEADERS_QTHREADS) $(PREFIX)/include/Qthreads

copy-openmp: mkdir
	mkdir -p $(PREFIX)/include/OpenMP
	cp $(COPY_FLAG) $(KOKKOS_HEADERS_OPENMP) $(PREFIX)/include/OpenMP

install: mkdir $(CONDITIONAL_COPIES) build-lib
	cp $(COPY_FLAG) $(NVCC_WRAPPER) $(PREFIX)/bin
	cp $(COPY_FLAG) $(KOKKOS_HEADERS_INCLUDE) $(PREFIX)/include
	cp $(COPY_FLAG) $(KOKKOS_HEADERS_INCLUDE_IMPL) $(PREFIX)/include/impl
	cp $(COPY_FLAG) Makefile.kokkos $(PREFIX)
	cp $(COPY_FLAG) kokkos.cmake $(PREFIX)
	cp $(COPY_FLAG) libkokkos.a $(PREFIX)/lib
	cp $(COPY_FLAG) KokkosCore_config.h $(PREFIX)/include

clean: kokkos-clean
	rm -f Makefile.kokkos
