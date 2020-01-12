ifndef KOKKOS_PATH
  MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
  KOKKOS_PATH = $(subst Makefile,,$(MAKEFILE_PATH))../..
endif

PREFIX ?= /usr/local/lib/kokkos

default: build-lib

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
  CXX ?= $(KOKKOS_PATH)/bin/nvcc_wrapper
else
  CXX ?= g++
endif

CXXFLAGS ?= -O3
LINK ?= $(CXX)
LDFLAGS ?=
CP = cp

include $(KOKKOS_PATH)/Makefile.kokkos
include $(KOKKOS_PATH)/core/src/Makefile.generate_header_lists
include $(KOKKOS_PATH)/core/src/Makefile.generate_build_files

CONDITIONAL_COPIES =

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
  CONDITIONAL_COPIES += copy-cuda
endif

ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)
  CONDITIONAL_COPIES += copy-threads
endif

ifeq ($(KOKKOS_INTERNAL_USE_QTHREADS), 1)
  CONDITIONAL_COPIES += copy-qthreads
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
  CONDITIONAL_COPIES += copy-openmp
endif

ifeq ($(KOKKOS_INTERNAL_USE_HPX), 1)
  CONDITIONAL_COPIES += copy-hpx
endif

ifeq ($(KOKKOS_INTERNAL_USE_ROCM), 1)
  CONDITIONAL_COPIES += copy-rocm
endif

ifeq ($(KOKKOS_OS),CYGWIN)
  COPY_FLAG = -u
endif
ifeq ($(KOKKOS_OS),Linux)
  COPY_FLAG = -u
endif
ifeq ($(KOKKOS_OS),Darwin)
  COPY_FLAG = 
  # If Homebrew coreutils is installed, its cp will have the -u option
  ifneq ("$(wildcard /usr/local/opt/coreutils/libexec/gnubin/cp)","")
    CP = /usr/local/opt/coreutils/libexec/gnubin/cp
    COPY_FLAG = -u
  endif
endif

ifeq ($(KOKKOS_DEBUG),"no")
  KOKKOS_DEBUG_CMAKE = OFF
else
  KOKKOS_DEBUG_CMAKE = ON
endif

build-lib: $(KOKKOS_LINK_DEPENDS)

mkdir:
	mkdir -p $(PREFIX)
	mkdir -p $(PREFIX)/bin
	mkdir -p $(PREFIX)/include
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/lib/pkgconfig
	mkdir -p $(PREFIX)/include/impl

copy-cuda: mkdir
	mkdir -p $(PREFIX)/include/Cuda
	$(CP) $(COPY_FLAG) $(KOKKOS_HEADERS_CUDA) $(PREFIX)/include/Cuda

copy-threads: mkdir
	mkdir -p $(PREFIX)/include/Threads
	$(CP) $(COPY_FLAG) $(KOKKOS_HEADERS_THREADS) $(PREFIX)/include/Threads

copy-qthreads: mkdir
	mkdir -p $(PREFIX)/include/Qthreads
	$(CP) $(COPY_FLAG) $(KOKKOS_HEADERS_QTHREADS) $(PREFIX)/include/Qthreads

copy-openmp: mkdir
	mkdir -p $(PREFIX)/include/OpenMP
	$(CP) $(COPY_FLAG) $(KOKKOS_HEADERS_OPENMP) $(PREFIX)/include/OpenMP

copy-hpx: mkdir
	mkdir -p $(PREFIX)/include/HPX
	$(CP) $(COPY_FLAG) $(KOKKOS_HEADERS_HPX) $(PREFIX)/include/HPX

copy-rocm: mkdir
	mkdir -p $(PREFIX)/include/ROCm
	$(CP) $(COPY_FLAG) $(KOKKOS_HEADERS_ROCM) $(PREFIX)/include/ROCm

install: mkdir $(CONDITIONAL_COPIES) build-lib generate_build_settings
	$(CP) $(COPY_FLAG) $(NVCC_WRAPPER) $(PREFIX)/bin
	$(CP) $(COPY_FLAG) $(KOKKOS_HEADERS_INCLUDE) $(PREFIX)/include
	$(CP) $(COPY_FLAG) $(KOKKOS_HEADERS_INCLUDE_IMPL) $(PREFIX)/include/impl
	$(CP) $(COPY_FLAG) $(KOKKOS_MAKEFILE)  $(PREFIX)
	$(CP) $(COPY_FLAG) $(KOKKOS_CMAKEFILE)  $(PREFIX)
	$(CP) $(COPY_FLAG) $(KOKKOS_PKGCONFIG)  $(PREFIX)/lib/pkgconfig
	$(CP) $(COPY_FLAG) libkokkos.a $(PREFIX)/lib
	$(CP) $(COPY_FLAG) $(KOKKOS_CONFIG_HEADER) $(PREFIX)/include

clean: kokkos-clean
	rm -f $(KOKKOS_MAKEFILE) $(KOKKOS_CMAKEFILE) $(KOKKOS_PKGCONFIG) 
