CUDA  = $(NVCC) $(CUDA_INCLUDE) $(CUDA_OPTS) -Icudpp_mini $(CUDA_ARCH) \
             $(CUDA_PRECISION)
CUDR  = $(CUDR_CPP) $(CUDR_OPTS) $(CUDA_PRECISION) $(CUDA_INCLUDE) \
         $(CUDPP_OPT)
CUDA_LINK = $(CUDA_LIB) -lcudart

GPU_LIB = $(LIB_DIR)/libgpu.a

# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
NVC_H  = $(wildcard ./geryon/nvc*.h) $(UCL_H)
NVD_H  = $(wildcard ./geryon/nvd*.h) $(UCL_H) lal_preprocessor.h
# Headers for Pair Stuff
PAIR_H  = lal_atom.h lal_answer.h lal_neighbor_shared.h \
          lal_neighbor.h lal_precision.h lal_device.h \
          lal_balance.h lal_pppm.h

ALL_H = $(NVD_H) $(PAIR_H)

EXECS = $(BIN_DIR)/nvc_get_devices
ifdef CUDPP_OPT
CUDPP = $(OBJ_DIR)/cudpp.o $(OBJ_DIR)/cudpp_plan.o \
        $(OBJ_DIR)/cudpp_maximal_launch.o $(OBJ_DIR)/cudpp_plan_manager.o \
        $(OBJ_DIR)/radixsort_app.cu_o $(OBJ_DIR)/scan_app.cu_o
endif
OBJS = $(OBJ_DIR)/lal_atom.o $(OBJ_DIR)/lal_ans.o \
       $(OBJ_DIR)/lal_neighbor.o $(OBJ_DIR)/lal_neighbor_shared.o \
       $(OBJ_DIR)/lal_device.o $(OBJ_DIR)/lal_base_atomic.o \
       $(OBJ_DIR)/lal_base_charge.o $(OBJ_DIR)/lal_base_ellipsoid.o \
       $(OBJ_DIR)/lal_pppm.o $(OBJ_DIR)/lal_pppm_ext.o \
       $(OBJ_DIR)/lal_gayberne.o $(OBJ_DIR)/lal_gayberne_ext.o \
       $(OBJ_DIR)/lal_re_squared.o $(OBJ_DIR)/lal_re_squared_ext.o \
       $(OBJ_DIR)/lal_lj.o $(OBJ_DIR)/lal_lj_ext.o \
       $(OBJ_DIR)/lal_lj96.o $(OBJ_DIR)/lal_lj96_ext.o \
       $(OBJ_DIR)/lal_lj_expand.o $(OBJ_DIR)/lal_lj_expand_ext.o \
       $(OBJ_DIR)/lal_lj_coul.o $(OBJ_DIR)/lal_lj_coul_ext.o \
       $(OBJ_DIR)/lal_lj_coul_long.o $(OBJ_DIR)/lal_lj_coul_long_ext.o \
       $(OBJ_DIR)/lal_lj_class2_long.o $(OBJ_DIR)/lal_lj_class2_long_ext.o \
       $(OBJ_DIR)/lal_coul_long.o $(OBJ_DIR)/lal_coul_long_ext.o \
       $(OBJ_DIR)/lal_morse.o $(OBJ_DIR)/lal_morse_ext.o \
       $(OBJ_DIR)/lal_charmm_long.o $(OBJ_DIR)/lal_charmm_long_ext.o \
       $(OBJ_DIR)/lal_cg_cmm.o $(OBJ_DIR)/lal_cg_cmm_ext.o \
       $(OBJ_DIR)/lal_cg_cmm_long.o $(OBJ_DIR)/lal_cg_cmm_long_ext.o \
       $(OBJ_DIR)/lal_eam.o $(OBJ_DIR)/lal_eam_ext.o \
       $(OBJ_DIR)/lal_buck.o $(OBJ_DIR)/lal_buck_ext.o \
       $(OBJ_DIR)/lal_buck_coul.o $(OBJ_DIR)/lal_buck_coul_ext.o \
       $(OBJ_DIR)/lal_buck_coul_long.o $(OBJ_DIR)/lal_buck_coul_long_ext.o \
       $(OBJ_DIR)/lal_table.o $(OBJ_DIR)/lal_table_ext.o \
       $(OBJ_DIR)/lal_yukawa.o $(OBJ_DIR)/lal_yukawa_ext.o
PTXS = $(OBJ_DIR)/device.ptx $(OBJ_DIR)/device_ptx.h \
       $(OBJ_DIR)/atom.ptx $(OBJ_DIR)/atom_ptx.h \
       $(OBJ_DIR)/neighbor_cpu.ptx $(OBJ_DIR)/neighbor_cpu_ptx.h \
       $(OBJ_DIR)/neighbor_gpu.ptx $(OBJ_DIR)/neighbor_gpu_ptx.h \
       $(OBJ_DIR)/pppm_f.ptx $(OBJ_DIR)/pppm_f_ptx.h \
       $(OBJ_DIR)/pppm_d.ptx $(OBJ_DIR)/pppm_d_ptx.h \
       $(OBJ_DIR)/ellipsoid_nbor.ptx $(OBJ_DIR)/ellipsoid_nbor_ptx.h \
       $(OBJ_DIR)/gayberne.ptx $(OBJ_DIR)/gayberne_lj.ptx \
       $(OBJ_DIR)/gayberne_ptx.h $(OBJ_DIR)/gayberne_lj_ptx.h \
       $(OBJ_DIR)/re_squared.ptx $(OBJ_DIR)/re_squared_lj.ptx \
       $(OBJ_DIR)/re_squared_ptx.h $(OBJ_DIR)/re_squared_lj_ptx.h \
       $(OBJ_DIR)/lj.ptx $(OBJ_DIR)/lj_ptx.h \
       $(OBJ_DIR)/lj96.ptx $(OBJ_DIR)/lj96_ptx.h \
       $(OBJ_DIR)/lj_expand.ptx $(OBJ_DIR)/lj_expand_ptx.h \
       $(OBJ_DIR)/lj_coul.ptx $(OBJ_DIR)/lj_coul_ptx.h \
       $(OBJ_DIR)/lj_coul_long.ptx $(OBJ_DIR)/lj_coul_long_ptx.h \
       $(OBJ_DIR)/lj_class2_long.ptx $(OBJ_DIR)/lj_class2_long_ptx.h \
       $(OBJ_DIR)/coul_long.ptx $(OBJ_DIR)/coul_long_ptx.h \
       $(OBJ_DIR)/morse.ptx $(OBJ_DIR)/morse_ptx.h \
       $(OBJ_DIR)/charmm_long.ptx $(OBJ_DIR)/charmm_long_ptx.h \
       $(OBJ_DIR)/cg_cmm.ptx $(OBJ_DIR)/cg_cmm_ptx.h \
       $(OBJ_DIR)/cg_cmm_long.ptx $(OBJ_DIR)/cg_cmm_long_ptx.h \
       $(OBJ_DIR)/eam.ptx $(OBJ_DIR)/eam_ptx.h \
       $(OBJ_DIR)/buck.ptx $(OBJ_DIR)/buck_ptx.h \
       $(OBJ_DIR)/buck_coul.ptx $(OBJ_DIR)/buck_coul_ptx.h \
       $(OBJ_DIR)/buck_coul_long.ptx $(OBJ_DIR)/buck_coul_long_ptx.h \
       $(OBJ_DIR)/table.ptx $(OBJ_DIR)/table_ptx.h \
       $(OBJ_DIR)/yukawa.ptx $(OBJ_DIR)/yukawa_ptx.h

all: $(GPU_LIB) $(EXECS)

$(OBJ_DIR)/cudpp.o: cudpp_mini/cudpp.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp.cpp -Icudpp_mini

$(OBJ_DIR)/cudpp_plan.o: cudpp_mini/cudpp_plan.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp_plan.cpp -Icudpp_mini

$(OBJ_DIR)/cudpp_maximal_launch.o: cudpp_mini/cudpp_maximal_launch.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp_maximal_launch.cpp -Icudpp_mini

$(OBJ_DIR)/cudpp_plan_manager.o: cudpp_mini/cudpp_plan_manager.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp_plan_manager.cpp -Icudpp_mini

$(OBJ_DIR)/radixsort_app.cu_o: cudpp_mini/radixsort_app.cu
	$(CUDA) -o $@ -c cudpp_mini/radixsort_app.cu

$(OBJ_DIR)/scan_app.cu_o: cudpp_mini/scan_app.cu
	$(CUDA) -o $@ -c cudpp_mini/scan_app.cu

$(OBJ_DIR)/atom.ptx: lal_atom.cu lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_atom.cu

$(OBJ_DIR)/atom_ptx.h: $(OBJ_DIR)/atom.ptx
	$(BSH) ./geryon/file_to_cstr.sh atom $(OBJ_DIR)/atom.ptx $(OBJ_DIR)/atom_ptx.h

$(OBJ_DIR)/lal_atom.o: lal_atom.cpp lal_atom.h $(NVD_H) $(OBJ_DIR)/atom_ptx.h
	$(CUDR) -o $@ -c lal_atom.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_ans.o: lal_answer.cpp lal_answer.h $(NVD_H)
	$(CUDR) -o $@ -c lal_answer.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/neighbor_cpu.ptx: lal_neighbor_cpu.cu lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_neighbor_cpu.cu

$(OBJ_DIR)/neighbor_cpu_ptx.h: $(OBJ_DIR)/neighbor_cpu.ptx
	$(BSH) ./geryon/file_to_cstr.sh neighbor_cpu $(OBJ_DIR)/neighbor_cpu.ptx $(OBJ_DIR)/neighbor_cpu_ptx.h

$(OBJ_DIR)/neighbor_gpu.ptx: lal_neighbor_gpu.cu lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_neighbor_gpu.cu

$(OBJ_DIR)/neighbor_gpu_ptx.h: $(OBJ_DIR)/neighbor_gpu.ptx
	$(BSH) ./geryon/file_to_cstr.sh neighbor_gpu $(OBJ_DIR)/neighbor_gpu.ptx $(OBJ_DIR)/neighbor_gpu_ptx.h

$(OBJ_DIR)/lal_neighbor_shared.o: lal_neighbor_shared.cpp lal_neighbor_shared.h $(OBJ_DIR)/neighbor_cpu_ptx.h $(OBJ_DIR)/neighbor_gpu_ptx.h $(NVD_H)
	$(CUDR) -o $@ -c lal_neighbor_shared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_neighbor.o: lal_neighbor.cpp lal_neighbor.h lal_neighbor_shared.h $(NVD_H)
	$(CUDR) -o $@ -c lal_neighbor.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/device.ptx: lal_device.cu lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_device.cu

$(OBJ_DIR)/device_ptx.h: $(OBJ_DIR)/device.ptx
	$(BSH) ./geryon/file_to_cstr.sh device $(OBJ_DIR)/device.ptx $(OBJ_DIR)/device_ptx.h

$(OBJ_DIR)/lal_device.o: lal_device.cpp lal_device.h $(ALL_H) $(OBJ_DIR)/device_ptx.h
	$(CUDR) -o $@ -c lal_device.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_base_atomic.o: $(ALL_H) lal_base_atomic.h lal_base_atomic.cpp
	$(CUDR) -o $@ -c lal_base_atomic.cpp

$(OBJ_DIR)/lal_base_charge.o: $(ALL_H) lal_base_charge.h lal_base_charge.cpp
	$(CUDR) -o $@ -c lal_base_charge.cpp

$(OBJ_DIR)/lal_base_ellipsoid.o: $(ALL_H) lal_base_ellipsoid.h lal_base_ellipsoid.cpp $(OBJ_DIR)/ellipsoid_nbor_ptx.h
	$(CUDR) -o $@ -c lal_base_ellipsoid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pppm_f.ptx: lal_pppm.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -Dgrdtyp=float -Dgrdtyp4=float4 -o $@ lal_pppm.cu

$(OBJ_DIR)/pppm_f_ptx.h: $(OBJ_DIR)/pppm_f.ptx
	$(BSH) ./geryon/file_to_cstr.sh pppm_f $(OBJ_DIR)/pppm_f.ptx $(OBJ_DIR)/pppm_f_ptx.h

$(OBJ_DIR)/pppm_d.ptx: lal_pppm.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -Dgrdtyp=double -Dgrdtyp4=double4 -o $@ lal_pppm.cu

$(OBJ_DIR)/pppm_d_ptx.h: $(OBJ_DIR)/pppm_d.ptx
	$(BSH) ./geryon/file_to_cstr.sh pppm_d $(OBJ_DIR)/pppm_d.ptx $(OBJ_DIR)/pppm_d_ptx.h

$(OBJ_DIR)/lal_pppm.o: $(ALL_H) lal_pppm.h lal_pppm.cpp $(OBJ_DIR)/pppm_f_ptx.h $(OBJ_DIR)/pppm_d_ptx.h
	$(CUDR) -o $@ -c lal_pppm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_pppm_ext.o: $(ALL_H) lal_pppm.h lal_pppm_ext.cpp
	$(CUDR) -o $@ -c lal_pppm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ellipsoid_nbor.ptx: lal_ellipsoid_nbor.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_ellipsoid_nbor.cu

$(OBJ_DIR)/ellipsoid_nbor_ptx.h: $(OBJ_DIR)/ellipsoid_nbor.ptx
	$(BSH) ./geryon/file_to_cstr.sh ellipsoid_nbor $(OBJ_DIR)/ellipsoid_nbor.ptx $(OBJ_DIR)/ellipsoid_nbor_ptx.h

$(OBJ_DIR)/gayberne.ptx: lal_gayberne.cu lal_precision.h lal_ellipsoid_extra.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_gayberne.cu

$(OBJ_DIR)/gayberne_lj.ptx: lal_gayberne_lj.cu lal_precision.h lal_ellipsoid_extra.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_gayberne_lj.cu

$(OBJ_DIR)/gayberne_ptx.h: $(OBJ_DIR)/gayberne.ptx
	$(BSH) ./geryon/file_to_cstr.sh gayberne $(OBJ_DIR)/gayberne.ptx $(OBJ_DIR)/gayberne_ptx.h

$(OBJ_DIR)/gayberne_lj_ptx.h: $(OBJ_DIR)/gayberne_lj.ptx
	$(BSH) ./geryon/file_to_cstr.sh gayberne_lj $(OBJ_DIR)/gayberne_lj.ptx $(OBJ_DIR)/gayberne_lj_ptx.h

$(OBJ_DIR)/lal_gayberne.o: $(ALL_H) lal_gayberne.h lal_gayberne.cpp $(OBJ_DIR)/gayberne_ptx.h $(OBJ_DIR)/gayberne_lj_ptx.h $(OBJ_DIR)/lal_base_ellipsoid.o
	$(CUDR) -o $@ -c lal_gayberne.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_gayberne_ext.o: $(ALL_H) $(OBJ_DIR)/lal_gayberne.o lal_gayberne_ext.cpp
	$(CUDR) -o $@ -c lal_gayberne_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/re_squared.ptx: lal_re_squared.cu lal_precision.h lal_ellipsoid_extra.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_re_squared.cu

$(OBJ_DIR)/re_squared_lj.ptx: lal_re_squared_lj.cu lal_precision.h lal_ellipsoid_extra.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_re_squared_lj.cu

$(OBJ_DIR)/re_squared_ptx.h: $(OBJ_DIR)/re_squared.ptx
	$(BSH) ./geryon/file_to_cstr.sh re_squared $(OBJ_DIR)/re_squared.ptx $(OBJ_DIR)/re_squared_ptx.h

$(OBJ_DIR)/re_squared_lj_ptx.h: $(OBJ_DIR)/re_squared_lj.ptx
	$(BSH) ./geryon/file_to_cstr.sh re_squared_lj $(OBJ_DIR)/re_squared_lj.ptx $(OBJ_DIR)/re_squared_lj_ptx.h

$(OBJ_DIR)/lal_re_squared.o: $(ALL_H) lal_re_squared.h lal_re_squared.cpp $(OBJ_DIR)/re_squared_ptx.h $(OBJ_DIR)/re_squared_lj_ptx.h $(OBJ_DIR)/lal_base_ellipsoid.o
	$(CUDR) -o $@ -c lal_re_squared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_re_squared_ext.o: $(ALL_H) $(OBJ_DIR)/lal_re_squared.o lal_re_squared_ext.cpp
	$(CUDR) -o $@ -c lal_re_squared_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj.ptx: lal_lj.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_lj.cu

$(OBJ_DIR)/lj_ptx.h: $(OBJ_DIR)/lj.ptx $(OBJ_DIR)/lj.ptx
	$(BSH) ./geryon/file_to_cstr.sh lj $(OBJ_DIR)/lj.ptx $(OBJ_DIR)/lj_ptx.h

$(OBJ_DIR)/lal_lj.o: $(ALL_H) lal_lj.h lal_lj.cpp $(OBJ_DIR)/lj_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_ext.o: $(ALL_H) lal_lj.h lal_lj_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_lj_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul.ptx: lal_lj_coul.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_lj_coul.cu

$(OBJ_DIR)/lj_coul_ptx.h: $(OBJ_DIR)/lj_coul.ptx $(OBJ_DIR)/lj_coul.ptx
	$(BSH) ./geryon/file_to_cstr.sh lj_coul $(OBJ_DIR)/lj_coul.ptx $(OBJ_DIR)/lj_coul_ptx.h

$(OBJ_DIR)/lal_lj_coul.o: $(ALL_H) lal_lj_coul.h lal_lj_coul.cpp $(OBJ_DIR)/lj_coul_ptx.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_ext.o: $(ALL_H) lal_lj_coul.h lal_lj_coul_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long.ptx: lal_lj_class2_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_lj_class2_long.cu

$(OBJ_DIR)/lj_class2_long_ptx.h: $(OBJ_DIR)/lj_class2_long.ptx $(OBJ_DIR)/lj_class2_long.ptx
	$(BSH) ./geryon/file_to_cstr.sh lj_class2_long $(OBJ_DIR)/lj_class2_long.ptx $(OBJ_DIR)/lj_class2_long_ptx.h

$(OBJ_DIR)/lal_lj_class2_long.o: $(ALL_H) lal_lj_class2_long.h lal_lj_class2_long.cpp $(OBJ_DIR)/lj_class2_long_ptx.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_class2_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_class2_long_ext.o: $(ALL_H) lal_lj_class2_long.h lal_lj_class2_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_class2_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_long.ptx: lal_coul_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_coul_long.cu

$(OBJ_DIR)/coul_long_ptx.h: $(OBJ_DIR)/coul_long.ptx $(OBJ_DIR)/coul_long.ptx
	$(BSH) ./geryon/file_to_cstr.sh coul_long $(OBJ_DIR)/coul_long.ptx $(OBJ_DIR)/coul_long_ptx.h

$(OBJ_DIR)/lal_coul_long.o: $(ALL_H) lal_coul_long.h lal_coul_long.cpp $(OBJ_DIR)/coul_long_ptx.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_long_ext.o: $(ALL_H) lal_coul_long.h lal_coul_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_long.ptx: lal_lj_coul_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_lj_coul_long.cu

$(OBJ_DIR)/lj_coul_long_ptx.h: $(OBJ_DIR)/lj_coul_long.ptx $(OBJ_DIR)/lj_coul_long.ptx
	$(BSH) ./geryon/file_to_cstr.sh lj_coul_long $(OBJ_DIR)/lj_coul_long.ptx $(OBJ_DIR)/lj_coul_long_ptx.h

$(OBJ_DIR)/lal_lj_coul_long.o: $(ALL_H) lal_lj_coul_long.h lal_lj_coul_long.cpp $(OBJ_DIR)/lj_coul_long_ptx.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_long_ext.o: $(ALL_H) lal_lj_coul_long.h lal_lj_coul_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse.ptx: lal_morse.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_morse.cu

$(OBJ_DIR)/morse_ptx.h: $(OBJ_DIR)/morse.ptx $(OBJ_DIR)/morse.ptx
	$(BSH) ./geryon/file_to_cstr.sh morse $(OBJ_DIR)/morse.ptx $(OBJ_DIR)/morse_ptx.h

$(OBJ_DIR)/lal_morse.o: $(ALL_H) lal_morse.h lal_morse.cpp $(OBJ_DIR)/morse_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_morse.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_morse_ext.o: $(ALL_H) lal_morse.h lal_morse_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_morse_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/charmm_long.ptx: lal_charmm_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_charmm_long.cu

$(OBJ_DIR)/charmm_long_ptx.h: $(OBJ_DIR)/charmm_long.ptx $(OBJ_DIR)/charmm_long.ptx
	$(BSH) ./geryon/file_to_cstr.sh charmm_long $(OBJ_DIR)/charmm_long.ptx $(OBJ_DIR)/charmm_long_ptx.h

$(OBJ_DIR)/lal_charmm_long.o: $(ALL_H) lal_charmm_long.h lal_charmm_long.cpp $(OBJ_DIR)/charmm_long_ptx.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_charmm_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_charmm_long_ext.o: $(ALL_H) lal_charmm_long.h lal_charmm_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_charmm_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96.ptx: lal_lj96.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_lj96.cu

$(OBJ_DIR)/lj96_ptx.h: $(OBJ_DIR)/lj96.ptx $(OBJ_DIR)/lj96.ptx
	$(BSH) ./geryon/file_to_cstr.sh lj96 $(OBJ_DIR)/lj96.ptx $(OBJ_DIR)/lj96_ptx.h

$(OBJ_DIR)/lal_lj96.o: $(ALL_H) lal_lj96.h lal_lj96.cpp $(OBJ_DIR)/lj96_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj96.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj96_ext.o: $(ALL_H) lal_lj96.h lal_lj96_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_lj96_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand.ptx: lal_lj_expand.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_lj_expand.cu

$(OBJ_DIR)/lj_expand_ptx.h: $(OBJ_DIR)/lj_expand.ptx $(OBJ_DIR)/lj_expand.ptx
	$(BSH) ./geryon/file_to_cstr.sh lj_expand $(OBJ_DIR)/lj_expand.ptx $(OBJ_DIR)/lj_expand_ptx.h

$(OBJ_DIR)/lal_lj_expand.o: $(ALL_H) lal_lj_expand.h lal_lj_expand.cpp $(OBJ_DIR)/lj_expand_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj_expand.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_expand_ext.o: $(ALL_H) lal_lj_expand.h lal_lj_expand_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_lj_expand_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm.ptx: lal_cg_cmm.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_cg_cmm.cu

$(OBJ_DIR)/cg_cmm_ptx.h: $(OBJ_DIR)/cg_cmm.ptx $(OBJ_DIR)/cg_cmm.ptx
	$(BSH) ./geryon/file_to_cstr.sh cg_cmm $(OBJ_DIR)/cg_cmm.ptx $(OBJ_DIR)/cg_cmm_ptx.h

$(OBJ_DIR)/lal_cg_cmm.o: $(ALL_H) lal_cg_cmm.h lal_cg_cmm.cpp $(OBJ_DIR)/cg_cmm_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_cg_cmm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_cg_cmm_ext.o: $(ALL_H) lal_cg_cmm.h lal_cg_cmm_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_cg_cmm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_long.ptx: lal_cg_cmm_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_cg_cmm_long.cu

$(OBJ_DIR)/cg_cmm_long_ptx.h: $(OBJ_DIR)/cg_cmm_long.ptx $(OBJ_DIR)/cg_cmm_long.ptx
	$(BSH) ./geryon/file_to_cstr.sh cg_cmm_long $(OBJ_DIR)/cg_cmm_long.ptx $(OBJ_DIR)/cg_cmm_long_ptx.h

$(OBJ_DIR)/lal_cg_cmm_long.o: $(ALL_H) lal_cg_cmm_long.h lal_cg_cmm_long.cpp $(OBJ_DIR)/cg_cmm_long_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_cg_cmm_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_cg_cmm_long_ext.o: $(ALL_H) lal_cg_cmm_long.h lal_cg_cmm_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_cg_cmm_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/eam.ptx: lal_eam.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_eam.cu
  
$(OBJ_DIR)/eam_ptx.h: $(OBJ_DIR)/eam.ptx $(OBJ_DIR)/eam.ptx
	$(BSH) ./geryon/file_to_cstr.sh eam $(OBJ_DIR)/eam.ptx $(OBJ_DIR)/eam_ptx.h
    
$(OBJ_DIR)/lal_eam.o: $(ALL_H) lal_eam.h lal_eam.cpp $(OBJ_DIR)/eam_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_eam.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_ext.o: $(ALL_H) lal_eam.h lal_eam_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_eam_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/buck.ptx: lal_buck.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_buck.cu
  
$(OBJ_DIR)/buck_ptx.h: $(OBJ_DIR)/buck.ptx $(OBJ_DIR)/buck.ptx
	$(BSH) ./geryon/file_to_cstr.sh buck $(OBJ_DIR)/buck.ptx $(OBJ_DIR)/buck_ptx.h
    
$(OBJ_DIR)/lal_buck.o: $(ALL_H) lal_buck.h lal_buck.cpp $(OBJ_DIR)/buck_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_buck.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_buck_ext.o: $(ALL_H) lal_buck.h lal_buck_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_buck_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/buck_coul.ptx: lal_buck_coul.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_buck_coul.cu
  
$(OBJ_DIR)/buck_coul_ptx.h: $(OBJ_DIR)/buck_coul.ptx $(OBJ_DIR)/buck_coul.ptx
	$(BSH) ./geryon/file_to_cstr.sh buck_coul $(OBJ_DIR)/buck_coul.ptx $(OBJ_DIR)/buck_coul_ptx.h
    
$(OBJ_DIR)/lal_buck_coul.o: $(ALL_H) lal_buck_coul.h lal_buck_coul.cpp $(OBJ_DIR)/buck_coul_ptx.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_buck_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_buck_coul_ext.o: $(ALL_H) lal_buck_coul.h lal_buck_coul_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_buck_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/buck_coul_long.ptx: lal_buck_coul_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_buck_coul_long.cu
  
$(OBJ_DIR)/buck_coul_long_ptx.h: $(OBJ_DIR)/buck_coul_long.ptx $(OBJ_DIR)/buck_coul_long.ptx
	$(BSH) ./geryon/file_to_cstr.sh buck_coul_long $(OBJ_DIR)/buck_coul_long.ptx $(OBJ_DIR)/buck_coul_long_ptx.h
    
$(OBJ_DIR)/lal_buck_coul_long.o: $(ALL_H) lal_buck_coul_long.h lal_buck_coul_long.cpp $(OBJ_DIR)/buck_coul_long_ptx.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_buck_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_buck_coul_long_ext.o: $(ALL_H) lal_buck_coul_long.h lal_buck_coul_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_buck_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/table.ptx: lal_table.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_table.cu
  
$(OBJ_DIR)/table_ptx.h: $(OBJ_DIR)/table.ptx $(OBJ_DIR)/table.ptx
	$(BSH) ./geryon/file_to_cstr.sh table $(OBJ_DIR)/table.ptx $(OBJ_DIR)/table_ptx.h
    
$(OBJ_DIR)/lal_table.o: $(ALL_H) lal_table.h lal_table.cpp $(OBJ_DIR)/table_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_table.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_table_ext.o: $(ALL_H) lal_table.h lal_table_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_table_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/yukawa.ptx: lal_yukawa.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lal_yukawa.cu
  
$(OBJ_DIR)/yukawa_ptx.h: $(OBJ_DIR)/yukawa.ptx $(OBJ_DIR)/yukawa.ptx
	$(BSH) ./geryon/file_to_cstr.sh yukawa $(OBJ_DIR)/yukawa.ptx $(OBJ_DIR)/yukawa_ptx.h
    
$(OBJ_DIR)/lal_yukawa.o: $(ALL_H) lal_yukawa.h lal_yukawa.cpp $(OBJ_DIR)/yukawa_ptx.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_yukawa.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_yukawa_ext.o: $(ALL_H) lal_yukawa.h lal_yukawa_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_yukawa_ext.cpp -I$(OBJ_DIR)

$(BIN_DIR)/nvc_get_devices: ./geryon/ucl_get_devices.cpp $(NVD_H)
	$(CUDR) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_CUDADR $(CUDA_LIB) -lcuda 

$(GPU_LIB): $(OBJS) $(CUDPP)
	$(AR) -crusv $(GPU_LIB) $(OBJS) $(CUDPP)

clean:
	rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(CUDPP) $(PTXS) *.linkinfo

veryclean: clean
	rm -rf *~ *.linkinfo

cleanlib:
	rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(PTXS) *.linkinfo
