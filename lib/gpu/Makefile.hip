# /* ----------------------------------------------------------------------   
#  Generic Linux Makefile for HIP
#     - export HIP_PLATFORM=hcc (or nvcc) before execution
#     - change HIP_ARCH for your GPU
# ------------------------------------------------------------------------- */

# this setting should match LAMMPS Makefile
# one of LAMMPS_SMALLBIG (default), LAMMPS_BIGBIG and LAMMPS_SMALLSMALL

LMP_INC = -DLAMMPS_SMALLBIG

# precision for GPU calculations
# -D_SINGLE_SINGLE  # Single precision for all calculations
# -D_DOUBLE_DOUBLE  # Double precision for all calculations
# -D_SINGLE_DOUBLE  # Accumulation of forces, etc. in double

HIP_PRECISION = -D_SINGLE_DOUBLE

HIP_OPTS = -O3 
HIP_HOST_OPTS = -Wno-deprecated-declarations
HIP_HOST_INCLUDE = 

# use device sort 
# requires linking with hipcc and hipCUB + (rocPRIM or CUB for AMD or Nvidia respectively)
HIP_HOST_OPTS += -DUSE_HIP_DEVICE_SORT 
# path to cub
HIP_HOST_INCLUDE += -I./
# path to hipcub
HIP_HOST_INCLUDE += -I$(HIP_PATH)/../include

# use mpi
HIP_HOST_OPTS += -DMPI_GERYON -DUCL_NO_EXIT
# this settings should match LAMMPS Makefile
MPI_COMP_OPTS = $(shell mpicxx --showme:compile)
MPI_LINK_OPTS = $(shell mpicxx --showme:link)
#MPI_COMP_OPTS += -I/usr/include/mpi  -DMPICH_IGNORE_CXX_SEEK -DOMPI_SKIP_MPICXX=1

HIP_PATH ?= $(wildcard /opt/rocm/hip)
HIP_PLATFORM=$(shell $(HIP_PATH)/bin/hipconfig --compiler)

ifeq (hcc,$(HIP_PLATFORM))
	HIP_OPTS  += -ffast-math
	# possible values: gfx803,gfx900,gfx906
	HIP_ARCH = gfx906
else ifeq (nvcc,$(HIP_PLATFORM))
	HIP_OPTS  += --use_fast_math
	HIP_ARCH = -gencode arch=compute_30,code=[sm_30,compute_30] -gencode arch=compute_32,code=[sm_32,compute_32] -gencode arch=compute_35,code=[sm_35,compute_35] \
		    -gencode arch=compute_50,code=[sm_50,compute_50] -gencode arch=compute_52,code=[sm_52,compute_52] -gencode arch=compute_53,code=[sm_53,compute_53]\
			-gencode arch=compute_60,code=[sm_60,compute_60] -gencode arch=compute_61,code=[sm_61,compute_61] -gencode arch=compute_62,code=[sm_62,compute_62]\
			-gencode arch=compute_70,code=[sm_70,compute_70] -gencode arch=compute_72,code=[sm_72,compute_72] -gencode arch=compute_75,code=[sm_75,compute_75] 
else
	$(error Specify HIP platform using 'export HIP_PLATFORM=(hcc,nvcc)')
endif

BIN_DIR = .
OBJ_DIR = ./obj
LIB_DIR = .
AR = ar
BSH = /bin/sh


# /* ----------------------------------------------------------------------   
#  				don't change section below without need			
# ------------------------------------------------------------------------- */

HIP_OPTS += -DUSE_HIP $(HIP_PRECISION)
HIP_GPU_OPTS += $(HIP_OPTS) -I./

ifeq (hcc,$(HIP_PLATFORM))
	HIP_HOST_OPTS += -fPIC
	HIP_GPU_CC  = $(HIP_PATH)/bin/hipcc --genco
	HIP_GPU_OPTS_S = -t="$(HIP_ARCH)" -f=\" 
	HIP_GPU_OPTS_E = \"
	HIP_KERNEL_SUFFIX = .cpp
	HIP_LIBS_TARGET = export HCC_AMDGPU_TARGET := $(HIP_ARCH)
	export HCC_AMDGPU_TARGET := $(HIP_ARCH)
else ifeq (nvcc,$(HIP_PLATFORM))
	HIP_GPU_CC  = $(HIP_PATH)/bin/hipcc --fatbin 
	HIP_GPU_OPTS += $(HIP_ARCH)
	HIP_GPU_SORT_ARCH = $(HIP_ARCH)
	# fix nvcc can't handle -pthread flag
	MPI_COMP_OPTS := $(subst -pthread,-Xcompiler -pthread,$(MPI_COMP_OPTS))
	MPI_LINK_OPTS := $(subst -pthread,-Xcompiler -pthread,$(MPI_LINK_OPTS))
endif

# hipcc is essential for device sort, because of hipcub is header only library and ROCm gpu code generation is deferred to the linking stage
HIP_HOST_CC = $(HIP_PATH)/bin/hipcc
HIP_HOST_OPTS += $(HIP_OPTS) $(MPI_COMP_OPTS) $(LMP_INC)
HIP_HOST_CC_CMD  = $(HIP_HOST_CC) $(HIP_HOST_OPTS) $(HIP_HOST_INCLUDE)

# sources

ALL_H  =  $(wildcard ./geryon/ucl*.h) $(wildcard ./geryon/hip*.h) $(wildcard ./lal_*.h)
SRCS := $(wildcard ./lal_*.cpp)
OBJS := $(subst ./,$(OBJ_DIR)/,$(SRCS:%.cpp=%.o))
CUS  := $(wildcard lal_*.cu)
CUHS := $(filter-out pppm_cubin.h, $(CUS:lal_%.cu=%_cubin.h)) pppm_f_cubin.h pppm_d_cubin.h
CUHS := $(addprefix $(OBJ_DIR)/, $(CUHS))

all: $(OBJ_DIR) $(CUHS) $(LIB_DIR)/libgpu.a $(BIN_DIR)/hip_get_devices

$(OBJ_DIR):
	mkdir -p $@    

# GPU kernels compilation 

$(OBJ_DIR)/pppm_f_cubin.h: lal_pppm.cu  $(ALL_H)
	@cp $< $(OBJ_DIR)/temp_pppm_f.cu$(HIP_KERNEL_SUFFIX)
	$(HIP_GPU_CC) $(HIP_GPU_OPTS_S) $(HIP_GPU_OPTS) -Dgrdtyp=float  -Dgrdtyp4=float4 $(HIP_GPU_OPTS_E)  -o $(OBJ_DIR)/pppm_f.cubin $(OBJ_DIR)/temp_pppm_f.cu$(HIP_KERNEL_SUFFIX)
	@xxd -i $(OBJ_DIR)/pppm_f.cubin $@
	@sed -i "s/[a-zA-Z0-9_]*pppm_f_cubin/pppm_f/g" $@
	@rm $(OBJ_DIR)/temp_pppm_f.cu$(HIP_KERNEL_SUFFIX) $(OBJ_DIR)/pppm_f.cubin

$(OBJ_DIR)/pppm_d_cubin.h: lal_pppm.cu  $(ALL_H)
	@cp $< $(OBJ_DIR)/temp_pppm_d.cu$(HIP_KERNEL_SUFFIX)
	$(HIP_GPU_CC) $(HIP_GPU_OPTS_S) $(HIP_GPU_OPTS) -Dgrdtyp=double -Dgrdtyp4=double4 $(HIP_GPU_OPTS_E)  -o $(OBJ_DIR)/pppm_d.cubin $(OBJ_DIR)/temp_pppm_d.cu$(HIP_KERNEL_SUFFIX)
	@xxd -i $(OBJ_DIR)/pppm_d.cubin $@
	@sed -i "s/[a-zA-Z0-9_]*pppm_d_cubin/pppm_d/g" $@
	@rm $(OBJ_DIR)/temp_pppm_d.cu$(HIP_KERNEL_SUFFIX) $(OBJ_DIR)/pppm_d.cubin

$(OBJ_DIR)/%_cubin.h: lal_%.cu  $(ALL_H)
	@cp $< $(OBJ_DIR)/temp_$*.cu$(HIP_KERNEL_SUFFIX)
	$(HIP_GPU_CC) $(HIP_GPU_OPTS_S) $(HIP_GPU_OPTS) $(HIP_GPU_OPTS_E)  -o $(OBJ_DIR)/$*.cubin $(OBJ_DIR)/temp_$*.cu$(HIP_KERNEL_SUFFIX)
	@xxd -i $(OBJ_DIR)/$*.cubin $@
	@sed -i "s/[a-zA-Z0-9_]*$*_cubin/$*/g" $@
	@rm $(OBJ_DIR)/temp_$*.cu$(HIP_KERNEL_SUFFIX) $(OBJ_DIR)/$*.cubin

# host sources compilation

$(OBJ_DIR)/lal_atom.o: lal_atom.cpp $(CUHS) $(ALL_H)
	$(HIP_HOST_CC_CMD) -o $@ -c $< -I$(OBJ_DIR) $(HIP_GPU_SORT_ARCH) 

$(OBJ_DIR)/lal_%.o: lal_%.cpp $(CUHS) $(ALL_H)
	$(HIP_HOST_CC_CMD) -o $@ -c $< -I$(OBJ_DIR)

# libgpu building

$(LIB_DIR)/libgpu.a: $(OBJS)
	$(AR) -crs $@ $(OBJS)
	echo -e "export HIP_PLATFORM := $(HIP_PLATFORM)\n$(HIP_LIBS_TARGET)" > 'Makefile.lammps'   

# test app building

$(BIN_DIR)/hip_get_devices: ./geryon/ucl_get_devices.cpp $(ALL_H)
	$(HIP_HOST_CC_CMD) -o $@ $< -DUCL_HIP $(MPI_LINK_OPTS)

clean:
	-rm -f $(BIN_DIR)/hip_get_devices $(LIB_DIR)/libgpu.a $(OBJS) $(OBJ_DIR)/temp_* $(CUHS)
