# makefile variables and settings related to configuring JIT with Lepton.

ENABLE_JIT=0
ifeq ($(shell uname -m),x86_64)
ENABLE_JIT=1
endif
ifeq ($(shell uname -m),amd64)
ENABLE_JIT=1
endif

LEPTON_INC = -I$(LEPTON_DIR)/include
LEPTON_DEF = -DLEPTON_BUILDING_STATIC_LIBRARY=1

ifeq ($(ENABLE_JIT),1)
LEPTON_INC += -I$(LEPTON_DIR)
LEPTON_DEF += -DLEPTON_USE_JIT=1 -DASMJIT_BUILD_X86=1 -DASMJIT_STATIC=1 -DASMJIT_BUILD_RELEASE=1
endif
