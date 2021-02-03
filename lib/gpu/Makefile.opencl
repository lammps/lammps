# /* ----------------------------------------------------------------------   
#  Generic Linux Makefile for OpenCL 
# ------------------------------------------------------------------------- */

# which file will be copied to Makefile.lammps

EXTRAMAKE = Makefile.lammps.opencl

# this setting should match LAMMPS Makefile
# one of LAMMPS_SMALLBIG (default), LAMMPS_BIGBIG and LAMMPS_SMALLSMALL

LMP_INC = -DLAMMPS_SMALLBIG

# precision for GPU calculations
# -D_SINGLE_SINGLE  # Single precision for all calculations
# -D_DOUBLE_DOUBLE  # Double precision for all calculations
# -D_SINGLE_DOUBLE  # Accumulation of forces, etc. in double

OCL_PREC = -D_SINGLE_DOUBLE

BIN_DIR = ./
OBJ_DIR = ./
LIB_DIR = ./
AR = ar
BSH = /bin/sh

# Compiler and linker settings

# OCL_TUNE = -DFERMI_OCL     # -- Uncomment for NVIDIA Fermi
# OCL_TUNE = -DKEPLER_OCL    # -- Uncomment for NVIDIA Kepler
# OCL_TUNE = -DCYPRESS_OCL   # -- Uncomment for AMD Cypress
OCL_TUNE = -DGENERIC_OCL   # -- Uncomment for generic device

OCL_INC = -I/usr/local/cuda/include  # Path to CL directory
OCL_CPP = mpic++ $(DEFAULT_DEVICE) -g -DMPI_GERYON -DUCL_NO_EXIT -DMPICH_IGNORE_CXX_SEEK $(LMP_INC) $(OCL_INC)
OCL_LINK = -lOpenCL
OCL  = $(OCL_CPP) $(OCL_PREC) $(OCL_TUNE) -DUSE_OPENCL

# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
OCL_H  = $(wildcard ./geryon/ocl*.h) $(UCL_H) lal_preprocessor.h
PRE1_H = lal_preprocessor.h lal_aux_fun1.h
ALL_H  =  $(OCL_H) $(wildcard ./lal_*.h)

# Source files
SRCS := $(wildcard ./lal_*.cpp)
OBJS := $(subst ./,$(OBJ_DIR)/,$(SRCS:%.cpp=%.o))
CUS  := $(wildcard lal_*.cu)
KERS := $(subst ./,$(OBJ_DIR)/,$(CUS:lal_%.cu=%_cl.h))
KERS := $(addprefix $(OBJ_DIR)/, $(KERS))

# targets

GPU_LIB = $(LIB_DIR)/libgpu.a

EXECS = $(BIN_DIR)/ocl_get_devices

all: $(OBJ_DIR) $(KERS) $(GPU_LIB) $(EXECS)

$(OBJ_DIR):
	mkdir -p $@

# device code compilation

$(OBJ_DIR)/%_cl.h: lal_%.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh $* $(PRE1_H) $< $@;

# host code compilation

$(OBJ_DIR)/lal_%.o: lal_%.cpp $(KERS)
	$(OCL) -o $@ -c $< -I$(OBJ_DIR)

# build libgpu.a

$(GPU_LIB): $(OBJS)
	$(AR) -crusv $(GPU_LIB) $(OBJS)
	@cp $(EXTRAMAKE) Makefile.lammps

# test app for querying device info

$(BIN_DIR)/ocl_get_devices: ./geryon/ucl_get_devices.cpp $(OCL_H)
	$(OCL) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_OPENCL $(OCL_LINK)

clean:
	-rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(KERS) *.linkinfo

veryclean: clean
	-rm -rf *~ *.linkinfo

cleanlib:
	-rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(KERS) *.linkinfo

