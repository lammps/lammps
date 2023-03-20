# Common headers for kernels
PRE1_H = lal_preprocessor.h lal_aux_fun1.h

# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
OCL_H  = $(wildcard ./geryon/ocl*.h) $(UCL_H) lal_precision.h

# Headers for Host files
HOST_H = lal_answer.h lal_atom.h lal_balance.h lal_base_atomic.h lal_base_amoeba.h \
         lal_base_charge.h lal_base_dipole.h lal_base_dpd.h \
         lal_base_ellipsoid.h lal_base_three.h lal_device.h lal_neighbor.h \
         lal_neighbor_shared.h lal_pre_ocl_config.h $(OCL_H)

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

# Compiler and linker

OCL  = $(OCL_CPP) $(OCL_PREC) $(OCL_TUNE) -DUSE_OPENCL

# device code compilation

$(OBJ_DIR)/atom_cl.h: lal_atom.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh atom lal_preprocessor.h lal_atom.cu $(OBJ_DIR)/atom_cl.h

$(OBJ_DIR)/neighbor_cpu_cl.h: lal_neighbor_cpu.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh neighbor_cpu lal_preprocessor.h lal_neighbor_cpu.cu $(OBJ_DIR)/neighbor_cpu_cl.h

$(OBJ_DIR)/neighbor_gpu_cl.h: lal_neighbor_gpu.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh neighbor_gpu lal_preprocessor.h lal_neighbor_gpu.cu $(OBJ_DIR)/neighbor_gpu_cl.h

$(OBJ_DIR)/device_cl.h: lal_device.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh device lal_preprocessor.h lal_device.cu $(OBJ_DIR)/device_cl.h

$(OBJ_DIR)/pppm_cl.h: lal_pppm.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh pppm lal_preprocessor.h lal_pppm.cu $(OBJ_DIR)/pppm_cl.h;

$(OBJ_DIR)/ellipsoid_nbor_cl.h: lal_ellipsoid_nbor.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh ellipsoid_nbor lal_preprocessor.h lal_ellipsoid_nbor.cu $(OBJ_DIR)/ellipsoid_nbor_cl.h

$(OBJ_DIR)/gayberne_cl.h: lal_gayberne.cu $(PRE1_H) lal_ellipsoid_extra.h
	$(BSH) ./geryon/file_to_cstr.sh gayberne $(PRE1_H) lal_ellipsoid_extra.h lal_gayberne.cu $(OBJ_DIR)/gayberne_cl.h;

$(OBJ_DIR)/gayberne_lj_cl.h: lal_gayberne_lj.cu $(PRE1_H) lal_ellipsoid_extra.h
	$(BSH) ./geryon/file_to_cstr.sh gayberne_lj $(PRE1_H) lal_ellipsoid_extra.h lal_gayberne_lj.cu $(OBJ_DIR)/gayberne_lj_cl.h;

$(OBJ_DIR)/re_squared_cl.h: lal_re_squared.cu $(PRE1_H) lal_ellipsoid_extra.h
	$(BSH) ./geryon/file_to_cstr.sh re_squared $(PRE1_H) lal_ellipsoid_extra.h lal_re_squared.cu $(OBJ_DIR)/re_squared_cl.h;

$(OBJ_DIR)/re_squared_lj_cl.h: lal_re_squared_lj.cu $(PRE1_H) lal_ellipsoid_extra.h
	$(BSH) ./geryon/file_to_cstr.sh re_squared_lj $(PRE1_H) lal_ellipsoid_extra.h lal_re_squared_lj.cu $(OBJ_DIR)/re_squared_lj_cl.h;

$(OBJ_DIR)/tersoff_cl.h: lal_tersoff.cu $(PRE1_H) lal_tersoff_extra.h
	$(BSH) ./geryon/file_to_cstr.sh tersoff $(PRE1_H) lal_tersoff_extra.h lal_tersoff.cu $(OBJ_DIR)/tersoff_cl.h;

$(OBJ_DIR)/tersoff_mod_cl.h: lal_tersoff_mod.cu $(PRE1_H) lal_tersoff_mod_extra.h
	$(BSH) ./geryon/file_to_cstr.sh tersoff_mod $(PRE1_H) lal_tersoff_mod_extra.h lal_tersoff_mod.cu $(OBJ_DIR)/tersoff_mod_cl.h;

$(OBJ_DIR)/tersoff_zbl_cl.h: lal_tersoff_zbl.cu $(PRE1_H) lal_tersoff_zbl_extra.h
	$(BSH) ./geryon/file_to_cstr.sh tersoff_zbl $(PRE1_H) lal_tersoff_zbl_extra.h lal_tersoff_zbl.cu $(OBJ_DIR)/tersoff_zbl_cl.h;

$(OBJ_DIR)/hippo_cl.h: lal_hippo.cu $(PRE1_H) lal_hippo_extra.h
	$(BSH) ./geryon/file_to_cstr.sh hippo $(PRE1_H) lal_hippo_extra.h lal_hippo.cu $(OBJ_DIR)/hippo_cl.h;

$(OBJ_DIR)/%_cl.h: lal_%.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh $* $(PRE1_H) $< $@;

# host code compilation

$(OBJ_DIR)/lal_answer.o: lal_answer.cpp $(HOST_H)
	$(OCL) -o $@ -c lal_answer.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dpd_tstat_ext.o: lal_dpd_tstat_ext.cpp lal_dpd.h $(HOST_H)
	$(OCL) -o $@ -c lal_dpd_tstat_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_alloy_ext.o: lal_eam_alloy_ext.cpp lal_eam.h $(HOST_H)
	$(OCL) -o $@ -c lal_eam_alloy_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_fs_ext.o: lal_eam_fs_ext.cpp lal_eam.h $(HOST_H)
	$(OCL) -o $@ -c lal_eam_fs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_neighbor.o: lal_neighbor.cpp $(HOST_H)
	$(OCL) -o $@ -c lal_neighbor.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_neighbor_shared.o: lal_neighbor_shared.cpp $(HOST_H)
	$(OCL) -o $@ -c lal_neighbor_shared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_%_ext.o: lal_%_ext.cpp lal_%.h $(HOST_H)
	$(OCL) -o $@ -c $< -I$(OBJ_DIR)

$(OBJ_DIR)/lal_base_%.o: lal_base_%.cpp $(HOST_H)
	$(OCL) -o $@ -c $< -I$(OBJ_DIR)

$(OBJ_DIR)/lal_%.o: lal_%.cpp %_cl.h $(HOST_H)
	$(OCL) -o $@ -c $< -I$(OBJ_DIR)

$(BIN_DIR)/ocl_get_devices: ./geryon/ucl_get_devices.cpp $(OCL_H)
	$(OCL) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_OPENCL $(OCL_LINK)

$(GPU_LIB): $(OBJS)
	$(AR) -crusv $(GPU_LIB) $(OBJS)
	@cp $(EXTRAMAKE) Makefile.lammps

clean:
	-rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(KERS) *.linkinfo

veryclean: clean
	-rm -rf *~ *.linkinfo

cleanlib:
	-rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(KERS) *.linkinfo

