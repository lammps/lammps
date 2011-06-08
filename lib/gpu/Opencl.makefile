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

OCL  = $(OCL_CPP) $(OCL_PREC) -DUSE_OPENCL
OCL_LIB = $(LIB_DIR)/libgpu.a
# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
OCL_H  = $(wildcard ./geryon/ocl*.h) $(UCL_H)
# Headers for Pair Stuff
PAIR_H  = pair_gpu_atom.h pair_gpu_ans.h pair_gpu_nbor_shared.h \
          pair_gpu_nbor.h pair_gpu_precision.h pair_gpu_device.h \
          pair_gpu_balance.h pppm_gpu_memory.h

ALL_H = $(OCL_H) $(PAIR_H)

EXECS = $(BIN_DIR)/ocl_get_devices
OBJS = $(OBJ_DIR)/pair_gpu_atom.o $(OBJ_DIR)/pair_gpu_ans.o \
       $(OBJ_DIR)/pair_gpu_nbor_shared.o $(OBJ_DIR)/pair_gpu_nbor.o \
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
       $(OBJ_DIR)/morse_gpu_memory.o $(OBJ_DIR)/morse_gpu.o \
       $(OBJ_DIR)/crml_gpu_memory.o $(OBJ_DIR)/crml_gpu.o \
       $(OBJ_DIR)/cmm_cut_gpu_memory.o $(OBJ_DIR)/cmm_cut_gpu.o \
       $(OBJ_DIR)/cmmc_long_gpu_memory.o $(OBJ_DIR)/cmmc_long_gpu.o 
KERS = $(OBJ_DIR)/pair_gpu_dev_cl.h $(OBJ_DIR)/pair_gpu_atom_cl.h \
       $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/pppm_gpu_cl.h \
       $(OBJ_DIR)/ellipsoid_nbor_cl.h $(OBJ_DIR)/gayberne_cl.h \
       $(OBJ_DIR)/re_squared_cl.h \
       $(OBJ_DIR)/lj_cut_gpu_cl.h $(OBJ_DIR)/lj96_cut_gpu_cl.h \
       $(OBJ_DIR)/lj_expand_gpu_cl.h $(OBJ_DIR)/ljc_cut_gpu_cl.h \
       $(OBJ_DIR)/ljcl_cut_gpu_cl.h $(OBJ_DIR)/lj_class2_long_cl.h \
       $(OBJ_DIR)/morse_gpu_cl.h \
       $(OBJ_DIR)/crml_gpu_cl.h $(OBJ_DIR)/cmm_cut_gpu_cl.h \
       $(OBJ_DIR)/cmmc_long_gpu_cl.h 

OCL_EXECS = $(BIN_DIR)/ocl_get_devices

all: $(OCL_LIB) $(EXECS)

$(OBJ_DIR)/pair_gpu_atom_cl.h: pair_gpu_atom_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh pair_gpu_atom_kernel.cu $(OBJ_DIR)/pair_gpu_atom_cl.h

$(OBJ_DIR)/pair_gpu_atom.o: pair_gpu_atom.cpp pair_gpu_atom.h $(OCL_H) $(OBJ_DIR)/pair_gpu_atom_cl.h
	$(OCL) -o $@ -c pair_gpu_atom.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_ans.o: pair_gpu_ans.cpp pair_gpu_ans.h $(OCL_H)
	$(OCL) -o $@ -c pair_gpu_ans.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_nbor_cl.h: pair_gpu_nbor_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh pair_gpu_nbor_kernel.cu $(OBJ_DIR)/pair_gpu_nbor_cl.h

$(OBJ_DIR)/pair_gpu_nbor_shared.o: pair_gpu_nbor_shared.cpp pair_gpu_nbor_shared.h $(OCL_H) $(OBJ_DIR)/pair_gpu_nbor_cl.h
	$(OCL) -o $@ -c pair_gpu_nbor_shared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_nbor.o: pair_gpu_nbor.cpp pair_gpu_nbor.h $(OCL_H) pair_gpu_nbor_shared.h
	$(OCL) -o $@ -c pair_gpu_nbor.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_dev_cl.h: pair_gpu_dev_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh pair_gpu_dev_kernel.cu $(OBJ_DIR)/pair_gpu_dev_cl.h

$(OBJ_DIR)/pair_gpu_device.o: pair_gpu_device.cpp pair_gpu_device.h $(ALL_H) $(OBJ_DIR)/pair_gpu_dev_cl.h
	$(OCL) -o $@ -c pair_gpu_device.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/atomic_gpu_memory.o: $(OCL_H) atomic_gpu_memory.h atomic_gpu_memory.cpp
	$(OCL) -o $@ -c atomic_gpu_memory.cpp

$(OBJ_DIR)/charge_gpu_memory.o: $(OCL_H) charge_gpu_memory.h charge_gpu_memory.cpp
	$(OCL) -o $@ -c charge_gpu_memory.cpp

$(OBJ_DIR)/base_ellipsoid.o: $(OCL_H) base_ellipsoid.h base_ellipsoid.cpp $(OBJ_DIR)/ellipsoid_nbor_cl.h
	$(OCL) -o $@ -c base_ellipsoid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pppm_gpu_cl.h: pppm_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh pppm_gpu_kernel.cu $(OBJ_DIR)/pppm_gpu_cl.h;

$(OBJ_DIR)/pppm_gpu_memory.o: $(ALL_H) pppm_gpu_memory.h pppm_gpu_memory.cpp  $(OBJ_DIR)/pppm_gpu_cl.h $(OBJ_DIR)/pppm_gpu_cl.h
	$(OCL) -o $@ -c pppm_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pppm_l_gpu.o: $(ALL_H) pppm_gpu_memory.h pppm_l_gpu.cpp
	$(OCL) -o $@ -c pppm_l_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ellipsoid_nbor_cl.h: ellipsoid_nbor.cu
	$(BSH) ./geryon/file_to_cstr.sh ellipsoid_nbor.cu $(OBJ_DIR)/ellipsoid_nbor_cl.h

$(OBJ_DIR)/gayberne_cl.h: gayberne.cu gayberne_lj.cu ellipsoid_extra.h
	cat ellipsoid_extra.h gayberne.cu > $(OBJ_DIR)/gayberne.tar; \
	cat ellipsoid_extra.h gayberne_lj.cu > $(OBJ_DIR)/gayberne_lj.tar; \
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/gayberne.tar $(OBJ_DIR)/gayberne_lj.tar $(OBJ_DIR)/gayberne_cl.h; \
	rm -f $(OBJ_DIR)/gayberne.tar $(OBJ_DIR)/gayberne_lj.tar

$(OBJ_DIR)/gayberne.o: $(ALL_H) gayberne.h gayberne.cpp $(OBJ_DIR)/gayberne_cl.h $(OBJ_DIR)/base_ellipsoid.o
	$(OCL) -o $@ -c gayberne.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/gayberne_ext.o: $(ALL_H) $(OBJ_DIR)/gayberne.o gayberne_ext.cpp
	$(OCL) -o $@ -c gayberne_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/re_squared_cl.h: re_squared.cu re_squared_lj.cu ellipsoid_extra.h
	cat ellipsoid_extra.h re_squared.cu > $(OBJ_DIR)/re_squared.tar; \
	cat ellipsoid_extra.h re_squared_lj.cu > $(OBJ_DIR)/re_squared_lj.tar; \
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/re_squared.tar $(OBJ_DIR)/re_squared_lj.tar $(OBJ_DIR)/re_squared_cl.h; \
	rm -f $(OBJ_DIR)/re_squared.tar $(OBJ_DIR)/re_squared_lj.tar

$(OBJ_DIR)/re_squared.o: $(ALL_H) re_squared.h re_squared.cpp $(OBJ_DIR)/re_squared_cl.h $(OBJ_DIR)/base_ellipsoid.o
	$(OCL) -o $@ -c re_squared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/re_squared_ext.o: $(ALL_H) $(OBJ_DIR)/re_squared.o re_squared_ext.cpp
	$(OCL) -o $@ -c re_squared_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cut_gpu_cl.h: lj_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_cut_gpu_kernel.cu $(OBJ_DIR)/lj_cut_gpu_cl.h;

$(OBJ_DIR)/lj_cut_gpu_memory.o: $(ALL_H) lj_cut_gpu_memory.h lj_cut_gpu_memory.cpp  $(OBJ_DIR)/lj_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/lj_cut_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c lj_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cut_gpu.o: $(ALL_H) lj_cut_gpu_memory.h lj_cut_gpu.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c lj_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljc_cut_gpu_cl.h: ljc_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh ljc_cut_gpu_kernel.cu $(OBJ_DIR)/ljc_cut_gpu_cl.h;

$(OBJ_DIR)/ljc_cut_gpu_memory.o: $(ALL_H) ljc_cut_gpu_memory.h ljc_cut_gpu_memory.cpp  $(OBJ_DIR)/ljc_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/ljc_cut_gpu_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c ljc_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljc_cut_gpu.o: $(ALL_H) ljc_cut_gpu_memory.h ljc_cut_gpu.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c ljc_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljcl_cut_gpu_cl.h: ljcl_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh ljcl_cut_gpu_kernel.cu $(OBJ_DIR)/ljcl_cut_gpu_cl.h;

$(OBJ_DIR)/ljcl_cut_gpu_memory.o: $(ALL_H) ljcl_cut_gpu_memory.h ljcl_cut_gpu_memory.cpp  $(OBJ_DIR)/ljcl_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c ljcl_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljcl_cut_gpu.o: $(ALL_H) ljcl_cut_gpu_memory.h ljcl_cut_gpu.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c ljcl_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long_cl.h: lj_class2_long.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_class2_long.cu $(OBJ_DIR)/lj_class2_long_cl.h;

$(OBJ_DIR)/lj_class2_long.o: $(ALL_H) lj_class2_long.h lj_class2_long.cpp  $(OBJ_DIR)/lj_class2_long_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c lj_class2_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long_ext.o: $(ALL_H) lj_class2_long.h lj_class2_long_ext.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c lj_class2_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse_gpu_cl.h: morse_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh morse_gpu_kernel.cu $(OBJ_DIR)/morse_gpu_cl.h;

$(OBJ_DIR)/morse_gpu_memory.o: $(ALL_H) morse_gpu_memory.h morse_gpu_memory.cpp  $(OBJ_DIR)/morse_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/morse_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c morse_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse_gpu.o: $(ALL_H) morse_gpu_memory.h morse_gpu.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c morse_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/crml_gpu_cl.h: crml_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh crml_gpu_kernel.cu $(OBJ_DIR)/crml_gpu_cl.h;

$(OBJ_DIR)/crml_gpu_memory.o: $(ALL_H) crml_gpu_memory.h crml_gpu_memory.cpp  $(OBJ_DIR)/crml_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/crml_gpu_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c crml_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/crml_gpu.o: $(ALL_H) crml_gpu_memory.h crml_gpu.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c crml_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_cut_gpu_cl.h: lj96_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh lj96_cut_gpu_kernel.cu $(OBJ_DIR)/lj96_cut_gpu_cl.h;

$(OBJ_DIR)/lj96_cut_gpu_memory.o: $(ALL_H) lj96_cut_gpu_memory.h lj96_cut_gpu_memory.cpp  $(OBJ_DIR)/lj96_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/lj96_cut_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c lj96_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_cut_gpu.o: $(ALL_H) lj96_cut_gpu_memory.h lj96_cut_gpu.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c lj96_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_gpu_cl.h: lj_expand_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_expand_gpu_kernel.cu $(OBJ_DIR)/lj_expand_gpu_cl.h;

$(OBJ_DIR)/lj_expand_gpu_memory.o: $(ALL_H) lj_expand_gpu_memory.h lj_expand_gpu_memory.cpp  $(OBJ_DIR)/lj_expand_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/lj_expand_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c lj_expand_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_gpu.o: $(ALL_H) lj_expand_gpu_memory.h lj_expand_gpu.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c lj_expand_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmm_cut_gpu_cl.h: cmm_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh cmm_cut_gpu_kernel.cu $(OBJ_DIR)/cmm_cut_gpu_cl.h;

$(OBJ_DIR)/cmm_cut_gpu_memory.o: $(ALL_H) cmm_cut_gpu_memory.h cmm_cut_gpu_memory.cpp  $(OBJ_DIR)/cmm_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/cmm_cut_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c cmm_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmm_cut_gpu.o: $(ALL_H) cmm_cut_gpu_memory.h cmm_cut_gpu.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c cmm_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmmc_long_gpu_cl.h: cmmc_long_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh cmmc_long_gpu_kernel.cu $(OBJ_DIR)/cmmc_long_gpu_cl.h;

$(OBJ_DIR)/cmmc_long_gpu_memory.o: $(ALL_H) cmmc_long_gpu_memory.h cmmc_long_gpu_memory.cpp  $(OBJ_DIR)/cmmc_long_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/cmmc_long_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c cmmc_long_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmmc_long_gpu.o: $(ALL_H) cmmc_long_gpu_memory.h cmmc_long_gpu.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c cmmc_long_gpu.cpp -I$(OBJ_DIR)

$(BIN_DIR)/ocl_get_devices: ./geryon/ucl_get_devices.cpp
	$(OCL) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_OPENCL $(OCL_LINK) 

$(OCL_LIB): $(OBJS) $(PTXS)
	$(AR) -crusv $(OCL_LIB) $(OBJS)

opencl: $(OCL_EXECS)

clean:
	rm -rf $(EXECS) $(OCL_EXECS) $(OCL_LIB) $(OBJS) $(KERS) *.linkinfo

veryclean: clean
	rm -rf *~ *.linkinfo

