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
          pair_gpu_balance.h pppm.h

ALL_H = $(OCL_H) $(PAIR_H)

EXECS = $(BIN_DIR)/ocl_get_devices
OBJS = $(OBJ_DIR)/pair_gpu_atom.o $(OBJ_DIR)/pair_gpu_ans.o \
       $(OBJ_DIR)/pair_gpu_nbor_shared.o $(OBJ_DIR)/pair_gpu_nbor.o \
       $(OBJ_DIR)/pair_gpu_device.o $(OBJ_DIR)/atomic_gpu_memory.o \
       $(OBJ_DIR)/charge_gpu_memory.o $(OBJ_DIR)/base_ellipsoid.o \
       $(OBJ_DIR)/pppm.o $(OBJ_DIR)/pppm_ext.o \
       $(OBJ_DIR)/gayberne.o $(OBJ_DIR)/gayberne_ext.o \
       $(OBJ_DIR)/re_squared.o $(OBJ_DIR)/re_squared_ext.o \
       $(OBJ_DIR)/lj.o $(OBJ_DIR)/lj_ext.o \
       $(OBJ_DIR)/lj96.o $(OBJ_DIR)/lj96_ext.o \
       $(OBJ_DIR)/lj_expand.o $(OBJ_DIR)/lj_expand_ext.o \
       $(OBJ_DIR)/lj_coul.o $(OBJ_DIR)/lj_coul_ext.o \
       $(OBJ_DIR)/lj_coul_long.o $(OBJ_DIR)/lj_coul_long_ext.o \
       $(OBJ_DIR)/lj_class2_long.o $(OBJ_DIR)/lj_class2_long_ext.o \
       $(OBJ_DIR)/morse.o $(OBJ_DIR)/morse_ext.o \
       $(OBJ_DIR)/charmm_long.o $(OBJ_DIR)/charmm_long_ext.o \
       $(OBJ_DIR)/cg_cmm.o $(OBJ_DIR)/cg_cmm_ext.o \
       $(OBJ_DIR)/cg_cmm_long.o $(OBJ_DIR)/cg_cmm_long_ext.o 
KERS = $(OBJ_DIR)/pair_gpu_dev_cl.h $(OBJ_DIR)/pair_gpu_atom_cl.h \
       $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/pppm_gpu_cl.h \
       $(OBJ_DIR)/ellipsoid_nbor_cl.h $(OBJ_DIR)/gayberne_cl.h \
       $(OBJ_DIR)/re_squared_cl.h \
       $(OBJ_DIR)/lj_ext_cl.h $(OBJ_DIR)/lj96_ext_cl.h \
       $(OBJ_DIR)/lj_expand_ext_cl.h $(OBJ_DIR)/lj_coul_ext_cl.h \
       $(OBJ_DIR)/lj_coul_long_ext_cl.h $(OBJ_DIR)/lj_class2_long_cl.h \
       $(OBJ_DIR)/morse_ext_cl.h \
       $(OBJ_DIR)/charmm_long_ext_cl.h $(OBJ_DIR)/cg_cmm_ext_cl.h \
       $(OBJ_DIR)/cg_cmm_long_ext_cl.h 

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

$(OBJ_DIR)/pppm_gpu_cl.h: pppm.cu
	$(BSH) ./geryon/file_to_cstr.sh pppm.cu $(OBJ_DIR)/pppm_gpu_cl.h;

$(OBJ_DIR)/pppm.o: $(ALL_H) pppm.h pppm.cpp  $(OBJ_DIR)/pppm_gpu_cl.h $(OBJ_DIR)/pppm_gpu_cl.h
	$(OCL) -o $@ -c pppm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pppm_ext.o: $(ALL_H) pppm.h pppm_ext.cpp
	$(OCL) -o $@ -c pppm_ext.cpp -I$(OBJ_DIR)

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

$(OBJ_DIR)/lj_ext_cl.h: lj.cu
	$(BSH) ./geryon/file_to_cstr.sh lj.cu $(OBJ_DIR)/lj_ext_cl.h;

$(OBJ_DIR)/lj.o: $(ALL_H) lj.h lj.cpp  $(OBJ_DIR)/lj_ext_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/lj_ext_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c lj.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_ext.o: $(ALL_H) lj.h lj_ext.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c lj_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_ext_cl.h: lj_coul.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_coul.cu $(OBJ_DIR)/lj_coul_ext_cl.h;

$(OBJ_DIR)/lj_coul.o: $(ALL_H) lj_coul.h lj_coul.cpp  $(OBJ_DIR)/lj_coul_ext_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/lj_coul_ext_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c lj_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_ext.o: $(ALL_H) lj_coul.h lj_coul_ext.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c lj_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_long_ext_cl.h: lj_coul_long.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_coul_long.cu $(OBJ_DIR)/lj_coul_long_ext_cl.h;

$(OBJ_DIR)/lj_coul_long.o: $(ALL_H) lj_coul_long.h lj_coul_long.cpp  $(OBJ_DIR)/lj_coul_long_ext_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c lj_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_long_ext.o: $(ALL_H) lj_coul_long.h lj_coul_long_ext.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c lj_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long_cl.h: lj_class2_long.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_class2_long.cu $(OBJ_DIR)/lj_class2_long_cl.h;

$(OBJ_DIR)/lj_class2_long.o: $(ALL_H) lj_class2_long.h lj_class2_long.cpp  $(OBJ_DIR)/lj_class2_long_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c lj_class2_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long_ext.o: $(ALL_H) lj_class2_long.h lj_class2_long_ext.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c lj_class2_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse_ext_cl.h: morse.cu
	$(BSH) ./geryon/file_to_cstr.sh morse.cu $(OBJ_DIR)/morse_ext_cl.h;

$(OBJ_DIR)/morse.o: $(ALL_H) morse.h morse.cpp  $(OBJ_DIR)/morse_ext_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/morse_ext_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c morse.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse_ext.o: $(ALL_H) morse.h morse_ext.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c morse_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/charmm_long_ext_cl.h: charmm_long.cu
	$(BSH) ./geryon/file_to_cstr.sh charmm_long.cu $(OBJ_DIR)/charmm_long_ext_cl.h;

$(OBJ_DIR)/charmm_long.o: $(ALL_H) charmm_long.h charmm_long.cpp  $(OBJ_DIR)/charmm_long_ext_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/charmm_long_ext_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c charmm_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/charmm_long_ext.o: $(ALL_H) charmm_long.h charmm_long_ext.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c charmm_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_ext_cl.h: lj96.cu
	$(BSH) ./geryon/file_to_cstr.sh lj96.cu $(OBJ_DIR)/lj96_ext_cl.h;

$(OBJ_DIR)/lj96.o: $(ALL_H) lj96.h lj96.cpp  $(OBJ_DIR)/lj96_ext_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/lj96_ext_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c lj96.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_ext.o: $(ALL_H) lj96.h lj96_ext.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c lj96_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_ext_cl.h: lj_expand.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_expand.cu $(OBJ_DIR)/lj_expand_ext_cl.h;

$(OBJ_DIR)/lj_expand.o: $(ALL_H) lj_expand.h lj_expand.cpp  $(OBJ_DIR)/lj_expand_ext_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/lj_expand_ext_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c lj_expand.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_ext.o: $(ALL_H) lj_expand.h lj_expand_ext.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c lj_expand_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_ext_cl.h: cg_cmm.cu
	$(BSH) ./geryon/file_to_cstr.sh cg_cmm.cu $(OBJ_DIR)/cg_cmm_ext_cl.h;

$(OBJ_DIR)/cg_cmm.o: $(ALL_H) cg_cmm.h cg_cmm.cpp  $(OBJ_DIR)/cg_cmm_ext_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/cg_cmm_ext_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c cg_cmm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_ext.o: $(ALL_H) cg_cmm.h cg_cmm_ext.cpp atomic_gpu_memory.h
	$(OCL) -o $@ -c cg_cmm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_long_ext_cl.h: cg_cmm_long.cu
	$(BSH) ./geryon/file_to_cstr.sh cg_cmm_long.cu $(OBJ_DIR)/cg_cmm_long_ext_cl.h;

$(OBJ_DIR)/cg_cmm_long.o: $(ALL_H) cg_cmm_long.h cg_cmm_long.cpp  $(OBJ_DIR)/cg_cmm_long_ext_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/cg_cmm_long_ext_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c cg_cmm_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_long_ext.o: $(ALL_H) cg_cmm_long.h cg_cmm_long_ext.cpp charge_gpu_memory.h
	$(OCL) -o $@ -c cg_cmm_long_ext.cpp -I$(OBJ_DIR)

$(BIN_DIR)/ocl_get_devices: ./geryon/ucl_get_devices.cpp
	$(OCL) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_OPENCL $(OCL_LINK) 

$(OCL_LIB): $(OBJS) $(PTXS)
	$(AR) -crusv $(OCL_LIB) $(OBJS)

opencl: $(OCL_EXECS)

clean:
	rm -rf $(EXECS) $(OCL_EXECS) $(OCL_LIB) $(OBJS) $(KERS) *.linkinfo

veryclean: clean
	rm -rf *~ *.linkinfo

