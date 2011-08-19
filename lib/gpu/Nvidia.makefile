# /* ----------------------------------------------------------------------   
#    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator       
#    http://lammps.sandia.gov, Sandia National Laboratories                   
#    Steve Plimpton, sjplimp@sandia.gov                                       
#                                                                             
#    Copyright (2003) Sandia Corporation.  Under the terms of Contract        
#    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains   
#    certain rights in this software.  This software is distributed under      
#    the GNU General Public License.                                          
#                                                                             
#    See the README file in the top-level LAMMPS directory.                   
# ------------------------------------------------------------------------- */
#                                                                             
# /* ----------------------------------------------------------------------   
#    Contributing authors: Mike Brown (ORNL), brownw@ornl.gov               
#                          Peng Wang (Nvidia), penwang@nvidia.com
#                          Inderaj Bains (NVIDIA), ibains@nvidia.com
#                          Paul Crozier (SNL), pscrozi@sandia.gov             
# ------------------------------------------------------------------------- */

CUDA  = $(NVCC) $(CUDA_INCLUDE) $(CUDA_OPTS) -Icudpp_mini $(CUDA_ARCH) \
             $(CUDA_PRECISION)
CUDR  = $(CUDR_CPP) $(CUDR_OPTS) $(CUDA_PRECISION) $(CUDA_INCLUDE) \
        -Icudpp_mini
CUDA_LINK = $(CUDA_LIB) -lcudart

GPU_LIB = $(LIB_DIR)/libgpu.a

# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
NVC_H  = $(wildcard ./geryon/nvc*.h) $(UCL_H)
NVD_H  = $(wildcard ./geryon/nvd*.h) $(UCL_H) nv_kernel_def.h
# Headers for Pair Stuff
PAIR_H  = pair_gpu_atom.h pair_gpu_ans.h pair_gpu_nbor_shared.h \
          pair_gpu_nbor.h pair_gpu_precision.h pair_gpu_device.h \
          pair_gpu_balance.h pppm_gpu_memory.h

ALL_H = $(NVD_H) $(PAIR_H)

EXECS = $(BIN_DIR)/nvc_get_devices
CUDPP = $(OBJ_DIR)/cudpp.o $(OBJ_DIR)/cudpp_plan.o \
        $(OBJ_DIR)/cudpp_maximal_launch.o $(OBJ_DIR)/cudpp_plan_manager.o \
        $(OBJ_DIR)/radixsort_app.cu_o $(OBJ_DIR)/scan_app.cu_o
OBJS = $(OBJ_DIR)/pair_gpu_atom.o $(OBJ_DIR)/pair_gpu_ans.o \
       $(OBJ_DIR)/pair_gpu_nbor.o $(OBJ_DIR)/pair_gpu_nbor_shared.o \
       $(OBJ_DIR)/pair_gpu_device.o $(OBJ_DIR)/atomic_gpu_memory.o \
       $(OBJ_DIR)/charge_gpu_memory.o $(OBJ_DIR)/base_ellipsoid.o \
       $(OBJ_DIR)/pppm_gpu_memory.o $(OBJ_DIR)/pppm_l_gpu.o \
       $(OBJ_DIR)/gayberne.o $(OBJ_DIR)/gayberne_ext.o \
       $(OBJ_DIR)/re_squared.o $(OBJ_DIR)/re_squared_ext.o \
       $(OBJ_DIR)/lj_cut_gpu_memory.o $(OBJ_DIR)/lj_cut_gpu.o \
       $(OBJ_DIR)/lj96_cut_gpu_memory.o $(OBJ_DIR)/lj96_cut_gpu.o \
       $(OBJ_DIR)/lj_expand_gpu_memory.o $(OBJ_DIR)/lj_expand_gpu.o \
       $(OBJ_DIR)/ljc_cut_gpu_memory.o $(OBJ_DIR)/ljc_cut_gpu.o \
       $(OBJ_DIR)/ljcl_cut_gpu_memory.o $(OBJ_DIR)/ljcl_cut_gpu.o \
       $(OBJ_DIR)/lj_class2_long.o $(OBJ_DIR)/lj_class2_long_ext.o \
       $(OBJ_DIR)/coul_long_gpu_memory.o $(OBJ_DIR)/coul_long_gpu.o \
       $(OBJ_DIR)/morse_gpu_memory.o $(OBJ_DIR)/morse_gpu.o \
       $(OBJ_DIR)/crml_gpu_memory.o $(OBJ_DIR)/crml_gpu.o \
       $(OBJ_DIR)/cmm_cut_gpu_memory.o $(OBJ_DIR)/cmm_cut_gpu.o \
       $(OBJ_DIR)/cmmc_long_gpu_memory.o $(OBJ_DIR)/cmmc_long_gpu.o \
       $(OBJ_DIR)/cmmc_msm_gpu_memory.o $(OBJ_DIR)/cmmc_msm_gpu.o \
       $(CUDPP)
PTXS = $(OBJ_DIR)/pair_gpu_dev_kernel.ptx \
       $(OBJ_DIR)/pair_gpu_atom_kernel.ptx $(OBJ_DIR)/pair_gpu_atom_ptx.h \
       $(OBJ_DIR)/pair_gpu_nbor_kernel.ptx $(OBJ_DIR)/pair_gpu_nbor_ptx.h \
       $(OBJ_DIR)/pair_gpu_build_kernel.ptx $(OBJ_DIR)/pair_gpu_build_ptx.h \
       $(OBJ_DIR)/pppm_f_gpu_kernel.ptx $(OBJ_DIR)/pppm_f_gpu_ptx.h \
       $(OBJ_DIR)/pppm_d_gpu_kernel.ptx $(OBJ_DIR)/pppm_d_gpu_ptx.h \
       $(OBJ_DIR)/ellipsoid_nbor.ptx $(OBJ_DIR)/ellipsoid_nbor_ptx.h \
       $(OBJ_DIR)/gayberne.ptx $(OBJ_DIR)/gayberne_lj.ptx \
       $(OBJ_DIR)/gayberne_ptx.h $(OBJ_DIR)/re_squared.ptx \
       $(OBJ_DIR)/re_squared_lj.ptx $(OBJ_DIR)/re_squared_ptx.h \
       $(OBJ_DIR)/lj_cut_gpu_kernel.ptx $(OBJ_DIR)/lj_cut_gpu_ptx.h \
       $(OBJ_DIR)/lj96_cut_gpu_kernel.ptx $(OBJ_DIR)/lj96_cut_gpu_ptx.h \
       $(OBJ_DIR)/lj_expand_gpu_kernel.ptx $(OBJ_DIR)/lj_expand_gpu_ptx.h \
       $(OBJ_DIR)/ljc_cut_gpu_kernel.ptx $(OBJ_DIR)/ljc_cut_gpu_ptx.h \
       $(OBJ_DIR)/ljcl_cut_gpu_kernel.ptx $(OBJ_DIR)/ljcl_cut_gpu_ptx.h \
       $(OBJ_DIR)/lj_class2_long.ptx $(OBJ_DIR)/lj_class2_long_ptx.h \
       $(OBJ_DIR)/coul_long_gpu_kernel.ptx $(OBJ_DIR)/coul_long_gpu_ptx.h \
       $(OBJ_DIR)/morse_gpu_kernel.ptx $(OBJ_DIR)/morse_gpu_ptx.h \
       $(OBJ_DIR)/crml_gpu_kernel.ptx $(OBJ_DIR)/crml_gpu_ptx.h \
       $(OBJ_DIR)/cmm_cut_gpu_kernel.ptx $(OBJ_DIR)/cmm_cut_gpu_ptx.h \
       $(OBJ_DIR)/cmmc_long_gpu_kernel.ptx $(OBJ_DIR)/cmmc_long_gpu_ptx.h \
       $(OBJ_DIR)/cmmc_msm_gpu_kernel.ptx $(OBJ_DIR)/cmmc_msm_gpu_ptx.h

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

$(OBJ_DIR)/pair_gpu_atom_kernel.ptx: pair_gpu_atom_kernel.cu
	$(CUDA) --ptx -DNV_KERNEL -o $@ pair_gpu_atom_kernel.cu

$(OBJ_DIR)/pair_gpu_atom_ptx.h: $(OBJ_DIR)/pair_gpu_atom_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/pair_gpu_atom_kernel.ptx $(OBJ_DIR)/pair_gpu_atom_ptx.h

$(OBJ_DIR)/pair_gpu_atom.o: pair_gpu_atom.cpp pair_gpu_atom.h $(NVD_H) $(OBJ_DIR)/pair_gpu_atom_ptx.h
	$(CUDR) -o $@ -c pair_gpu_atom.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_ans.o: pair_gpu_ans.cpp pair_gpu_ans.h $(NVD_H)
	$(CUDR) -o $@ -c pair_gpu_ans.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_nbor_kernel.ptx: pair_gpu_nbor_kernel.cu
	$(CUDA) --ptx -DNV_KERNEL -o $@ pair_gpu_nbor_kernel.cu

$(OBJ_DIR)/pair_gpu_nbor_ptx.h: $(OBJ_DIR)/pair_gpu_nbor_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/pair_gpu_nbor_kernel.ptx $(OBJ_DIR)/pair_gpu_nbor_ptx.h

$(OBJ_DIR)/pair_gpu_build_kernel.ptx: pair_gpu_build_kernel.cu
	$(CUDA) --ptx -DNV_KERNEL -o $@ pair_gpu_build_kernel.cu

$(OBJ_DIR)/pair_gpu_build_ptx.h: $(OBJ_DIR)/pair_gpu_build_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/pair_gpu_build_kernel.ptx $(OBJ_DIR)/pair_gpu_build_ptx.h

$(OBJ_DIR)/pair_gpu_nbor_shared.o: pair_gpu_nbor_shared.cpp pair_gpu_nbor_shared.h $(OBJ_DIR)/pair_gpu_nbor_ptx.h $(OBJ_DIR)/pair_gpu_build_ptx.h $(NVD_H)
	$(CUDR) -o $@ -c pair_gpu_nbor_shared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_nbor.o: pair_gpu_nbor.cpp pair_gpu_nbor.h pair_gpu_nbor_shared.h $(NVD_H)
	$(CUDR) -o $@ -c pair_gpu_nbor.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_dev_kernel.ptx: pair_gpu_dev_kernel.cu
	$(CUDA) --ptx -DNV_KERNEL -o $@ pair_gpu_dev_kernel.cu

$(OBJ_DIR)/pair_gpu_dev_ptx.h: $(OBJ_DIR)/pair_gpu_dev_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/pair_gpu_dev_kernel.ptx $(OBJ_DIR)/pair_gpu_dev_ptx.h

$(OBJ_DIR)/pair_gpu_device.o: pair_gpu_device.cpp pair_gpu_device.h $(ALL_H) $(OBJ_DIR)/pair_gpu_dev_ptx.h
	$(CUDR) -o $@ -c pair_gpu_device.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/atomic_gpu_memory.o: $(ALL_H) atomic_gpu_memory.h atomic_gpu_memory.cpp
	$(CUDR) -o $@ -c atomic_gpu_memory.cpp

$(OBJ_DIR)/charge_gpu_memory.o: $(ALL_H) charge_gpu_memory.h charge_gpu_memory.cpp
	$(CUDR) -o $@ -c charge_gpu_memory.cpp

$(OBJ_DIR)/base_ellipsoid.o: $(ALL_H) base_ellipsoid.h base_ellipsoid.cpp $(OBJ_DIR)/ellipsoid_nbor_ptx.h
	$(CUDR) -o $@ -c base_ellipsoid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pppm_f_gpu_kernel.ptx: pppm_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -Dgrdtyp=float -Dgrdtyp4=float4 -o $@ pppm_gpu_kernel.cu

$(OBJ_DIR)/pppm_f_gpu_ptx.h: $(OBJ_DIR)/pppm_f_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/pppm_f_gpu_kernel.ptx $(OBJ_DIR)/pppm_f_gpu_ptx.h

$(OBJ_DIR)/pppm_d_gpu_kernel.ptx: pppm_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -Dgrdtyp=double -Dgrdtyp4=double4 -o $@ pppm_gpu_kernel.cu

$(OBJ_DIR)/pppm_d_gpu_ptx.h: $(OBJ_DIR)/pppm_d_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/pppm_d_gpu_kernel.ptx $(OBJ_DIR)/pppm_d_gpu_ptx.h

$(OBJ_DIR)/pppm_gpu_memory.o: $(ALL_H) pppm_gpu_memory.h pppm_gpu_memory.cpp $(OBJ_DIR)/pppm_f_gpu_ptx.h $(OBJ_DIR)/pppm_d_gpu_ptx.h
	$(CUDR) -o $@ -c pppm_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pppm_l_gpu.o: $(ALL_H) pppm_gpu_memory.h pppm_l_gpu.cpp
	$(CUDR) -o $@ -c pppm_l_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ellipsoid_nbor.ptx: ellipsoid_nbor.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ ellipsoid_nbor.cu

$(OBJ_DIR)/ellipsoid_nbor_ptx.h: $(OBJ_DIR)/ellipsoid_nbor.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/ellipsoid_nbor.ptx $(OBJ_DIR)/ellipsoid_nbor_ptx.h

$(OBJ_DIR)/gayberne.ptx: gayberne.cu pair_gpu_precision.h ellipsoid_extra.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ gayberne.cu

$(OBJ_DIR)/gayberne_lj.ptx: gayberne_lj.cu pair_gpu_precision.h ellipsoid_extra.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ gayberne_lj.cu

$(OBJ_DIR)/gayberne_ptx.h: $(OBJ_DIR)/gayberne.ptx $(OBJ_DIR)/gayberne_lj.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/gayberne.ptx $(OBJ_DIR)/gayberne_lj.ptx $(OBJ_DIR)/gayberne_ptx.h

$(OBJ_DIR)/gayberne.o: $(ALL_H) gayberne.h gayberne.cpp $(OBJ_DIR)/gayberne_ptx.h $(OBJ_DIR)/base_ellipsoid.o
	$(CUDR) -o $@ -c gayberne.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/gayberne_ext.o: $(ALL_H) $(OBJ_DIR)/gayberne.o gayberne_ext.cpp
	$(CUDR) -o $@ -c gayberne_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/re_squared.ptx: re_squared.cu pair_gpu_precision.h ellipsoid_extra.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ re_squared.cu

$(OBJ_DIR)/re_squared_lj.ptx: re_squared_lj.cu pair_gpu_precision.h ellipsoid_extra.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ re_squared_lj.cu

$(OBJ_DIR)/re_squared_ptx.h: $(OBJ_DIR)/re_squared.ptx $(OBJ_DIR)/re_squared_lj.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/re_squared.ptx $(OBJ_DIR)/re_squared_lj.ptx $(OBJ_DIR)/re_squared_ptx.h

$(OBJ_DIR)/re_squared.o: $(ALL_H) re_squared.h re_squared.cpp $(OBJ_DIR)/re_squared_ptx.h $(OBJ_DIR)/base_ellipsoid.o
	$(CUDR) -o $@ -c re_squared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/re_squared_ext.o: $(ALL_H) $(OBJ_DIR)/re_squared.o re_squared_ext.cpp
	$(CUDR) -o $@ -c re_squared_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cut_gpu_kernel.ptx: lj_cut_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lj_cut_gpu_kernel.cu

$(OBJ_DIR)/lj_cut_gpu_ptx.h: $(OBJ_DIR)/lj_cut_gpu_kernel.ptx $(OBJ_DIR)/lj_cut_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/lj_cut_gpu_kernel.ptx $(OBJ_DIR)/lj_cut_gpu_ptx.h

$(OBJ_DIR)/lj_cut_gpu_memory.o: $(ALL_H) lj_cut_gpu_memory.h lj_cut_gpu_memory.cpp $(OBJ_DIR)/lj_cut_gpu_ptx.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(CUDR) -o $@ -c lj_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cut_gpu.o: $(ALL_H) lj_cut_gpu_memory.h lj_cut_gpu.cpp atomic_gpu_memory.h
	$(CUDR) -o $@ -c lj_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljc_cut_gpu_kernel.ptx: ljc_cut_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ ljc_cut_gpu_kernel.cu

$(OBJ_DIR)/ljc_cut_gpu_ptx.h: $(OBJ_DIR)/ljc_cut_gpu_kernel.ptx $(OBJ_DIR)/ljc_cut_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/ljc_cut_gpu_kernel.ptx $(OBJ_DIR)/ljc_cut_gpu_ptx.h

$(OBJ_DIR)/ljc_cut_gpu_memory.o: $(ALL_H) ljc_cut_gpu_memory.h ljc_cut_gpu_memory.cpp $(OBJ_DIR)/ljc_cut_gpu_ptx.h $(OBJ_DIR)/charge_gpu_memory.o
	$(CUDR) -o $@ -c ljc_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljc_cut_gpu.o: $(ALL_H) ljc_cut_gpu_memory.h ljc_cut_gpu.cpp charge_gpu_memory.h
	$(CUDR) -o $@ -c ljc_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long.ptx: lj_class2_long.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lj_class2_long.cu

$(OBJ_DIR)/lj_class2_long_ptx.h: $(OBJ_DIR)/lj_class2_long.ptx $(OBJ_DIR)/lj_class2_long.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/lj_class2_long.ptx $(OBJ_DIR)/lj_class2_long_ptx.h

$(OBJ_DIR)/lj_class2_long.o: $(ALL_H) lj_class2_long.h lj_class2_long.cpp $(OBJ_DIR)/lj_class2_long_ptx.h $(OBJ_DIR)/charge_gpu_memory.o
	$(CUDR) -o $@ -c lj_class2_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long_ext.o: $(ALL_H) lj_class2_long.h lj_class2_long_ext.cpp charge_gpu_memory.h
	$(CUDR) -o $@ -c lj_class2_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljcl_cut_gpu_kernel.ptx: ljcl_cut_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ ljcl_cut_gpu_kernel.cu

$(OBJ_DIR)/ljcl_cut_gpu_ptx.h: $(OBJ_DIR)/ljcl_cut_gpu_kernel.ptx $(OBJ_DIR)/ljcl_cut_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/ljcl_cut_gpu_kernel.ptx $(OBJ_DIR)/ljcl_cut_gpu_ptx.h

$(OBJ_DIR)/ljcl_cut_gpu_memory.o: $(ALL_H) ljcl_cut_gpu_memory.h ljcl_cut_gpu_memory.cpp $(OBJ_DIR)/ljcl_cut_gpu_ptx.h $(OBJ_DIR)/charge_gpu_memory.o
	$(CUDR) -o $@ -c ljcl_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljcl_cut_gpu.o: $(ALL_H) ljcl_cut_gpu_memory.h ljcl_cut_gpu.cpp charge_gpu_memory.h
	$(CUDR) -o $@ -c ljcl_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_long_gpu_kernel.ptx: coul_long_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ coul_long_gpu_kernel.cu

$(OBJ_DIR)/coul_long_gpu_ptx.h: $(OBJ_DIR)/coul_long_gpu_kernel.ptx $(OBJ_DIR)/coul_long_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/coul_long_gpu_kernel.ptx $(OBJ_DIR)/coul_long_gpu_ptx.h

$(OBJ_DIR)/coul_long_gpu_memory.o: $(ALL_H) coul_long_gpu_memory.h coul_long_gpu_memory.cpp $(OBJ_DIR)/coul_long_gpu_ptx.h $(OBJ_DIR)/charge_gpu_memory.o
	$(CUDR) -o $@ -c coul_long_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_long_gpu.o: $(ALL_H) coul_long_gpu_memory.h coul_long_gpu.cpp charge_gpu_memory.h
	$(CUDR) -o $@ -c coul_long_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse_gpu_kernel.ptx: morse_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ morse_gpu_kernel.cu

$(OBJ_DIR)/morse_gpu_ptx.h: $(OBJ_DIR)/morse_gpu_kernel.ptx $(OBJ_DIR)/morse_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/morse_gpu_kernel.ptx $(OBJ_DIR)/morse_gpu_ptx.h

$(OBJ_DIR)/morse_gpu_memory.o: $(ALL_H) morse_gpu_memory.h morse_gpu_memory.cpp $(OBJ_DIR)/morse_gpu_ptx.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(CUDR) -o $@ -c morse_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse_gpu.o: $(ALL_H) morse_gpu_memory.h morse_gpu.cpp atomic_gpu_memory.h
	$(CUDR) -o $@ -c morse_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/crml_gpu_kernel.ptx: crml_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ crml_gpu_kernel.cu

$(OBJ_DIR)/crml_gpu_ptx.h: $(OBJ_DIR)/crml_gpu_kernel.ptx $(OBJ_DIR)/crml_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/crml_gpu_kernel.ptx $(OBJ_DIR)/crml_gpu_ptx.h

$(OBJ_DIR)/crml_gpu_memory.o: $(ALL_H) crml_gpu_memory.h crml_gpu_memory.cpp $(OBJ_DIR)/crml_gpu_ptx.h $(OBJ_DIR)/charge_gpu_memory.o
	$(CUDR) -o $@ -c crml_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/crml_gpu.o: $(ALL_H) crml_gpu_memory.h crml_gpu.cpp charge_gpu_memory.h
	$(CUDR) -o $@ -c crml_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_cut_gpu_kernel.ptx: lj96_cut_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lj96_cut_gpu_kernel.cu

$(OBJ_DIR)/lj96_cut_gpu_ptx.h: $(OBJ_DIR)/lj96_cut_gpu_kernel.ptx $(OBJ_DIR)/lj96_cut_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/lj96_cut_gpu_kernel.ptx $(OBJ_DIR)/lj96_cut_gpu_ptx.h

$(OBJ_DIR)/lj96_cut_gpu_memory.o: $(ALL_H) lj96_cut_gpu_memory.h lj96_cut_gpu_memory.cpp $(OBJ_DIR)/lj96_cut_gpu_ptx.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(CUDR) -o $@ -c lj96_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_cut_gpu.o: $(ALL_H) lj96_cut_gpu_memory.h lj96_cut_gpu.cpp atomic_gpu_memory.h
	$(CUDR) -o $@ -c lj96_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_gpu_kernel.ptx: lj_expand_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ lj_expand_gpu_kernel.cu

$(OBJ_DIR)/lj_expand_gpu_ptx.h: $(OBJ_DIR)/lj_expand_gpu_kernel.ptx $(OBJ_DIR)/lj_expand_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/lj_expand_gpu_kernel.ptx $(OBJ_DIR)/lj_expand_gpu_ptx.h

$(OBJ_DIR)/lj_expand_gpu_memory.o: $(ALL_H) lj_expand_gpu_memory.h lj_expand_gpu_memory.cpp $(OBJ_DIR)/lj_expand_gpu_ptx.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(CUDR) -o $@ -c lj_expand_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_gpu.o: $(ALL_H) lj_expand_gpu_memory.h lj_expand_gpu.cpp atomic_gpu_memory.h
	$(CUDR) -o $@ -c lj_expand_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmm_cut_gpu_kernel.ptx: cmm_cut_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ cmm_cut_gpu_kernel.cu

$(OBJ_DIR)/cmm_cut_gpu_ptx.h: $(OBJ_DIR)/cmm_cut_gpu_kernel.ptx $(OBJ_DIR)/cmm_cut_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/cmm_cut_gpu_kernel.ptx $(OBJ_DIR)/cmm_cut_gpu_ptx.h

$(OBJ_DIR)/cmm_cut_gpu_memory.o: $(ALL_H) cmm_cut_gpu_memory.h cmm_cut_gpu_memory.cpp $(OBJ_DIR)/cmm_cut_gpu_ptx.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(CUDR) -o $@ -c cmm_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmm_cut_gpu.o: $(ALL_H) cmm_cut_gpu_memory.h cmm_cut_gpu.cpp atomic_gpu_memory.h
	$(CUDR) -o $@ -c cmm_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmmc_long_gpu_kernel.ptx: cmmc_long_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ cmmc_long_gpu_kernel.cu

$(OBJ_DIR)/cmmc_long_gpu_ptx.h: $(OBJ_DIR)/cmmc_long_gpu_kernel.ptx $(OBJ_DIR)/cmmc_long_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/cmmc_long_gpu_kernel.ptx $(OBJ_DIR)/cmmc_long_gpu_ptx.h

$(OBJ_DIR)/cmmc_long_gpu_memory.o: $(ALL_H) cmmc_long_gpu_memory.h cmmc_long_gpu_memory.cpp $(OBJ_DIR)/cmmc_long_gpu_ptx.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(CUDR) -o $@ -c cmmc_long_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmmc_long_gpu.o: $(ALL_H) cmmc_long_gpu_memory.h cmmc_long_gpu.cpp charge_gpu_memory.h
	$(CUDR) -o $@ -c cmmc_long_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmmc_msm_gpu_kernel.ptx: cmmc_msm_gpu_kernel.cu pair_gpu_precision.h
	$(CUDA) --ptx -DNV_KERNEL -o $@ cmmc_msm_gpu_kernel.cu

$(OBJ_DIR)/cmmc_msm_gpu_ptx.h: $(OBJ_DIR)/cmmc_msm_gpu_kernel.ptx $(OBJ_DIR)/cmmc_msm_gpu_kernel.ptx
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/cmmc_msm_gpu_kernel.ptx $(OBJ_DIR)/cmmc_msm_gpu_ptx.h

$(OBJ_DIR)/cmmc_msm_gpu_memory.o: $(ALL_H) cmmc_msm_gpu_memory.h cmmc_msm_gpu_memory.cpp $(OBJ_DIR)/cmmc_msm_gpu_ptx.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(CUDR) -o $@ -c cmmc_msm_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmmc_msm_gpu.o: $(ALL_H) cmmc_msm_gpu_memory.h cmmc_msm_gpu.cpp charge_gpu_memory.h
	$(CUDR) -o $@ -c cmmc_msm_gpu.cpp -I$(OBJ_DIR)

$(BIN_DIR)/nvc_get_devices: ./geryon/ucl_get_devices.cpp $(NVC_H)
	$(CUDR) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_CUDART $(CUDA_LINK) 

$(GPU_LIB): $(OBJS)
	$(AR) -crusv $(GPU_LIB) $(OBJS)

clean:
	rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(PTXS) *.linkinfo

veryclean: clean
	rm -rf *~ *.linkinfo
