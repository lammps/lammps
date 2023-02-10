# Common headers for kernels
PRE1_H = lal_preprocessor.h lal_aux_fun1.h

# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
NVD_H  = $(wildcard ./geryon/nvd*.h) $(UCL_H) lal_preprocessor.h \
         lal_pre_cuda_hip.h

# Headers for Host files
HOST_H = lal_answer.h lal_atom.h lal_balance.h lal_base_atomic.h lal_base_amoeba.h \
         lal_base_charge.h lal_base_dipole.h lal_base_dpd.h \
         lal_base_ellipsoid.h lal_base_three.h lal_device.h lal_neighbor.h \
         lal_neighbor_shared.h lal_pre_ocl_config.h $(NVD_H)
         
# Source files
SRCS := $(wildcard ./lal_*.cpp)
OBJS := $(subst ./,$(OBJ_DIR)/,$(SRCS:%.cpp=%.o))
CUS  := $(wildcard lal_*.cu)
CUHS := $(filter-out pppm_cubin.h, $(CUS:lal_%.cu=%_cubin.h)) pppm_f_cubin.h pppm_d_cubin.h
CUHS := $(addprefix $(OBJ_DIR)/, $(CUHS))

ifdef CUDPP_OPT
CUDPP = $(OBJ_DIR)/cudpp.o $(OBJ_DIR)/cudpp_plan.o \
        $(OBJ_DIR)/cudpp_maximal_launch.o $(OBJ_DIR)/cudpp_plan_manager.o \
        $(OBJ_DIR)/radixsort_app.cu_o $(OBJ_DIR)/scan_app.cu_o
endif

# targets

GPU_LIB = $(LIB_DIR)/libgpu.a

EXECS = $(BIN_DIR)/nvc_get_devices

all: $(OBJ_DIR) $(CUHS) $(GPU_LIB) $(EXECS)

$(OBJ_DIR):
	mkdir -p $@

# Compilers and linkers

CUDA  = $(NVCC) $(CUDA_INCLUDE) $(CUDA_OPTS) $(CUDA_ARCH) \
             $(CUDA_PRECISION)
CUDR  = $(CUDR_CPP) $(CUDR_OPTS) $(CUDA_PRECISION) $(CUDA_INCLUDE) \
         $(CUDPP_OPT)
CUDA_LINK = $(CUDA_LIB) -lcudart

BIN2C = $(CUDA_HOME)/bin/bin2c

# device code compilation

$(OBJ_DIR)/pppm_f.cubin: lal_pppm.cu lal_precision.h lal_preprocessor.h \
                         lal_pre_cuda_hip.h
	$(CUDA) --fatbin -DNV_KERNEL -Dgrdtyp=float -Dgrdtyp4=float4 -o $@ lal_pppm.cu

$(OBJ_DIR)/pppm_f_cubin.h: $(OBJ_DIR)/pppm_f.cubin
	$(BIN2C) -c -n pppm_f $(OBJ_DIR)/pppm_f.cubin > $(OBJ_DIR)/pppm_f_cubin.h

$(OBJ_DIR)/pppm_d.cubin: lal_pppm.cu lal_precision.h lal_preprocessor.h \
                         lal_pre_cuda_hip.h
	$(CUDA) --fatbin -DNV_KERNEL -Dgrdtyp=double -Dgrdtyp4=double4 -o $@ lal_pppm.cu

$(OBJ_DIR)/pppm_d_cubin.h: $(OBJ_DIR)/pppm_d.cubin
	$(BIN2C) -c -n pppm_d $(OBJ_DIR)/pppm_d.cubin > $(OBJ_DIR)/pppm_d_cubin.h

$(OBJ_DIR)/%_cubin.h: lal_%.cu $(PRE1_H)
	$(CUDA) --fatbin -DNV_KERNEL -o $(OBJ_DIR)/$*.cubin $(OBJ_DIR)/lal_$*.cu
	$(BIN2C) -c -n $* $(OBJ_DIR)/$*.cubin > $@

# host code compilation

$(OBJ_DIR)/lal_answer.o: lal_answer.cpp $(HOST_H)
	$(CUDR) -o $@ -c lal_answer.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dpd_tstat_ext.o: lal_dpd_tstat_ext.cpp lal_dpd.h $(HOST_H)
	$(CUDR) -o $@ -c lal_dpd_tstat_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_alloy_ext.o: lal_eam_alloy_ext.cpp lal_eam.h $(HOST_H)
	$(CUDR) -o $@ -c lal_eam_alloy_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_fs_ext.o: lal_eam_fs_ext.cpp lal_eam.h $(HOST_H)
	$(CUDR) -o $@ -c lal_eam_fs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_neighbor.o: lal_neighbor.cpp $(HOST_H)
	$(CUDR) -o $@ -c lal_neighbor.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_neighbor_shared.o: lal_neighbor_shared.cpp $(HOST_H)
	$(CUDR) -o $@ -c lal_neighbor_shared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_pppm.o: lal_pppm.cpp pppm_f_cubin.h pppm_d_cubin.h $(HOST_H)
	$(CUDR) -o $@ -c lal_pppm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_%_ext.o: lal_%_ext.cpp lal_%.h $(HOST_H)
	$(CUDR) -o $@ -c $< -I$(OBJ_DIR)

$(OBJ_DIR)/lal_base_%.o: lal_base_%.cpp $(HOST_H)
	$(CUDR) -o $@ -c $< -I$(OBJ_DIR)

$(OBJ_DIR)/lal_%.o: lal_%.cpp %_cubin.h $(HOST_H)
	$(CUDR) -o $@ -c $< -I$(OBJ_DIR)

#ifdef CUDPP_OPT
$(OBJ_DIR)/cudpp.o: cudpp_mini/cudpp.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp.cpp -Icudpp_mini

$(OBJ_DIR)/cudpp_plan.o: cudpp_mini/cudpp_plan.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp_plan.cpp -Icudpp_mini

$(OBJ_DIR)/cudpp_maximal_launch.o: cudpp_mini/cudpp_maximal_launch.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp_maximal_launch.cpp -Icudpp_mini

$(OBJ_DIR)/cudpp_plan_manager.o: cudpp_mini/cudpp_plan_manager.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp_plan_manager.cpp -Icudpp_mini

$(OBJ_DIR)/radixsort_app.cu_o: cudpp_mini/radixsort_app.cu
	$(CUDA) -o $@ -c cudpp_mini/radixsort_app.cu -Icudpp_mini

$(OBJ_DIR)/scan_app.cu_o: cudpp_mini/scan_app.cu
	$(CUDA) -o $@ -c cudpp_mini/scan_app.cu -Icudpp_mini
#endif

# build libgpu.a

$(GPU_LIB): $(OBJS) $(CUDPP)
	$(AR) -crusv $(GPU_LIB) $(OBJS) $(CUDPP)
	@cp $(EXTRAMAKE) Makefile.lammps

# test app for querying device info

$(BIN_DIR)/nvc_get_devices: ./geryon/ucl_get_devices.cpp $(NVD_H)
	$(CUDR) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_CUDADR $(CUDA_LIB) -lcuda 

clean:
	-rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(CUDPP) $(CUHS) *.cubin *.linkinfo

veryclean: clean
	-rm -rf *~ *.linkinfo

cleanlib:
	-rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(CUHS) *.linkinfo
