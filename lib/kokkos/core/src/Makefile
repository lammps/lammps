KOKKOS_PATH = ../..

PREFIX ?= /usr/local/lib/kokkos

default: messages build-lib
	echo "End Build"


include $(KOKKOS_PATH)/Makefile.kokkos

ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
	CXX = $(NVCC_WRAPPER)
	CXXFLAGS ?= -O3
	LINK = $(NVCC_WRAPPER)
	LINKFLAGS ?= 
else
	CXX ?= g++
	CXXFLAGS ?= -O3
	LINK ?= g++
	LINKFLAGS ?=  
endif

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

ifeq ($(KOKKOS_INTERNAL_USE_QTHREAD), 1)
	KOKKOS_HEADERS_QTHREAD += $(wildcard $(KOKKOS_PATH)/core/src/Qthread/*.hpp)
	CONDITIONAL_COPIES += copy-qthread
endif

ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
	KOKKOS_HEADERS_OPENMP += $(wildcard $(KOKKOS_PATH)/core/src/OpenMP/*.hpp)
	CONDITIONAL_COPIES += copy-openmp
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
	echo "KOKKOS_CPP_DEPENDS = $(KOKKOS_CPP_DEPENDS)" >> Makefile.kokkos
	echo "KOKKOS_CXXFLAGS = $(KOKKOS_CXXFLAGS)" >> Makefile.kokkos
	echo "KOKKOS_CPPFLAGS = $(KOKKOS_CPPFLAGS)" >> Makefile.kokkos
	echo "KOKKOS_LINK_DEPENDS  = $(KOKKOS_LINK_DEPENDS)" >> Makefile.kokkos
	echo "KOKKOS_LIBS = $(KOKKOS_LIBS)" >> Makefile.kokkos
	echo "KOKKOS_LDFLAGS = $(KOKKOS_LDFLAGS)" >> Makefile.kokkos
	sed \
		-e 's|$(KOKKOS_PATH)/core/src|$(PREFIX)/include|g' \
		-e 's|$(KOKKOS_PATH)/containers/src|$(PREFIX)/include|g' \
		-e 's|$(KOKKOS_PATH)/algorithms/src|$(PREFIX)/include|g' \
		-e 's|-L$(PWD)|-L$(PREFIX)/lib|g' \
		-e 's|= libkokkos.a|= $(PREFIX)/lib/libkokkos.a|g' \
		-e 's|= KokkosCore_config.h|= $(PREFIX)/include/KokkosCore_config.h|g' Makefile.kokkos \
		> Makefile.kokkos.tmp
	mv -f Makefile.kokkos.tmp Makefile.kokkos

build-lib: build-makefile-kokkos $(KOKKOS_LINK_DEPENDS)

mkdir: 
	mkdir -p $(PREFIX)
	mkdir -p $(PREFIX)/bin
	mkdir -p $(PREFIX)/include
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include/impl

copy-cuda: mkdir
	mkdir -p $(PREFIX)/include/Cuda
	cp $(KOKKOS_HEADERS_CUDA) $(PREFIX)/include/Cuda

copy-threads: mkdir
	mkdir -p $(PREFIX)/include/Threads
	cp $(KOKKOS_HEADERS_THREADS) $(PREFIX)/include/Threads

copy-qthread: mkdir
	mkdir -p $(PREFIX)/include/Qthread
	cp $(KOKKOS_HEADERS_QTHREAD) $(PREFIX)/include/Qthread

copy-openmp: mkdir
	mkdir -p $(PREFIX)/include/OpenMP
	cp $(KOKKOS_HEADERS_OPENMP) $(PREFIX)/include/OpenMP

install: mkdir $(CONDITIONAL_COPIES) build-lib 
	cp $(NVCC_WRAPPER) $(PREFIX)/bin
	cp $(KOKKOS_HEADERS_INCLUDE) $(PREFIX)/include
	cp $(KOKKOS_HEADERS_INCLUDE_IMPL) $(PREFIX)/include/impl
	cp Makefile.kokkos $(PREFIX)
	cp libkokkos.a $(PREFIX)/lib
	cp KokkosCore_config.h $(PREFIX)/include



clean: kokkos-clean
	rm Makefile.kokkos




