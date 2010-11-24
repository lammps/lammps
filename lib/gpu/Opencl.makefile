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
#                          Paul Crozier (SNL), pscrozi@sandia.gov             
# ------------------------------------------------------------------------- */

OCL  = $(OCL_CPP) $(OCL_PREC) -DUSE_OPENCL
OCL_LIB = $(LIB_DIR)/libgpu.a
# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
OCL_H  = $(wildcard ./geryon/ocl*.h) $(UCL_H)
# Headers for Pair Stuff
PAIR_H  = pair_gpu_atom.h pair_gpu_nbor.h pair_gpu_precision.h \
          pair_gpu_device.h pair_gpu_balance.h

ALL_H = $(OCL_H) $(PAIR_H)

EXECS = $(BIN_DIR)/ocl_get_devices
OBJS = $(OBJ_DIR)/pair_gpu_atom.o $(OBJ_DIR)/pair_gpu_nbor.o \
       $(OBJ_DIR)/pair_gpu_device.o $(OBJ_DIR)/atomic_gpu_memory.o \
       $(OBJ_DIR)/charge_gpu_memory.o \
       $(OBJ_DIR)/gb_gpu_memory.o $(OBJ_DIR)/gb_gpu.o \
       $(OBJ_DIR)/lj_cut_gpu_memory.o $(OBJ_DIR)/lj_cut_gpu.o \
       $(OBJ_DIR)/lj96_cut_gpu_memory.o $(OBJ_DIR)/lj96_cut_gpu.o \
       $(OBJ_DIR)/ljc_cut_gpu_memory.o $(OBJ_DIR)/ljc_cut_gpu.o \
       $(OBJ_DIR)/ljcl_cut_gpu_memory.o $(OBJ_DIR)/ljcl_cut_gpu.o \
       $(OBJ_DIR)/cmm_cut_gpu_memory.o $(OBJ_DIR)/cmm_cut_gpu.o \
       $(OBJ_DIR)/cmmc_long_gpu_memory.o $(OBJ_DIR)/cmmc_long_gpu.o 
KERS = $(OBJ_DIR)/pair_gpu_atom_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h \
       $(OBJ_DIR)/gb_gpu_nbor_cl.h $(OBJ_DIR)/gb_gpu_cl.h \
       $(OBJ_DIR)/lj_cut_gpu_cl.h $(OBJ_DIR)/lj96_cut_gpu_cl.h \
       $(OBJ_DIR)/ljc_cut_gpu_cl.h $(OBJ_DIR)/ljcl_cut_gpu_cl.h \
       $(OBJ_DIR)/cmm_cut_gpu_cl.h $(OBJ_DIR)/cmmc_long_gpu_cl.h 
       
OCL_EXECS = $(BIN_DIR)/ocl_get_devices

all: $(OCL_LIB) $(EXECS)

$(OBJ_DIR)/pair_gpu_atom_cl.h: pair_gpu_atom_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh pair_gpu_atom_kernel.cu $(OBJ_DIR)/pair_gpu_atom_cl.h

$(OBJ_DIR)/pair_gpu_atom.o: pair_gpu_atom.cpp pair_gpu_atom.h $(OCL_H) $(OBJ_DIR)/pair_gpu_atom_cl.h
	$(OCL) -o $@ -c pair_gpu_atom.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_nbor_cl.h: pair_gpu_nbor_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh pair_gpu_nbor_kernel.cu $(OBJ_DIR)/pair_gpu_nbor_cl.h

$(OBJ_DIR)/pair_gpu_nbor.o: pair_gpu_nbor.cpp pair_gpu_nbor.h $(OCL_H) $(OBJ_DIR)/pair_gpu_nbor_cl.h
	$(OCL) -o $@ -c pair_gpu_nbor.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pair_gpu_device.o: pair_gpu_device.cpp pair_gpu_device.h $(OCL_H)
	$(OCL) -o $@ -c pair_gpu_device.cpp

$(OBJ_DIR)/atomic_gpu_memory.o: $(OCL_H) atomic_gpu_memory.h atomic_gpu_memory.cpp
	$(OCL) -o $@ -c atomic_gpu_memory.cpp

$(OBJ_DIR)/charge_gpu_memory.o: $(OCL_H) charge_gpu_memory.h charge_gpu_memory.cpp
	$(OCL) -o $@ -c charge_gpu_memory.cpp

$(OBJ_DIR)/gb_gpu_nbor_cl.h: gb_gpu_kernel_nbor.cu
	$(BSH) ./geryon/file_to_cstr.sh gb_gpu_kernel_nbor.cu $(OBJ_DIR)/gb_gpu_nbor_cl.h

$(OBJ_DIR)/gb_gpu_cl.h: gb_gpu_kernel.cu gb_gpu_kernel_lj.cu gb_gpu_extra.h
	cat gb_gpu_extra.h gb_gpu_kernel.cu > $(OBJ_DIR)/gb_gpu_kernel.tar; \
	cat gb_gpu_extra.h gb_gpu_kernel_lj.cu > $(OBJ_DIR)/gb_gpu_kernel_lj.tar; \
	$(BSH) ./geryon/file_to_cstr.sh $(OBJ_DIR)/gb_gpu_kernel.tar $(OBJ_DIR)/gb_gpu_kernel_lj.tar $(OBJ_DIR)/gb_gpu_cl.h; \
	rm -f $(OBJ_DIR)/gb_gpu_kernel.tar $(OBJ_DIR)/gb_gpu_kernel_lj.tar

$(OBJ_DIR)/gb_gpu_memory.o: $(ALL_H) gb_gpu_memory.h gb_gpu_memory.cpp $(OBJ_DIR)/gb_gpu_nbor_cl.h $(OBJ_DIR)/gb_gpu_cl.h
	$(OCL) -o $@ -c gb_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/gb_gpu.o: $(ALL_H) gb_gpu_memory.h gb_gpu.cpp
	$(OCL) -o $@ -c gb_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cut_gpu_cl.h: lj_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_cut_gpu_kernel.cu $(OBJ_DIR)/lj_cut_gpu_cl.h;

$(OBJ_DIR)/lj_cut_gpu_memory.o: $(ALL_H) lj_cut_gpu_memory.h lj_cut_gpu_memory.cpp  $(OBJ_DIR)/lj_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/lj_cut_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c lj_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cut_gpu.o: $(ALL_H) lj_cut_gpu_memory.h lj_cut_gpu.cpp
	$(OCL) -o $@ -c lj_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljc_cut_gpu_cl.h: ljc_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh ljc_cut_gpu_kernel.cu $(OBJ_DIR)/ljc_cut_gpu_cl.h;

$(OBJ_DIR)/ljc_cut_gpu_memory.o: $(ALL_H) ljc_cut_gpu_memory.h ljc_cut_gpu_memory.cpp  $(OBJ_DIR)/ljc_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/ljc_cut_gpu_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c ljc_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljc_cut_gpu.o: $(ALL_H) ljc_cut_gpu_memory.h ljc_cut_gpu.cpp
	$(OCL) -o $@ -c ljc_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljcl_cut_gpu_cl.h: ljcl_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh ljcl_cut_gpu_kernel.cu $(OBJ_DIR)/ljcl_cut_gpu_cl.h;

$(OBJ_DIR)/ljcl_cut_gpu_memory.o: $(ALL_H) ljcl_cut_gpu_memory.h ljcl_cut_gpu_memory.cpp  $(OBJ_DIR)/ljcl_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/ljcl_cut_gpu_cl.h $(OBJ_DIR)/charge_gpu_memory.o
	$(OCL) -o $@ -c ljcl_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ljcl_cut_gpu.o: $(ALL_H) ljcl_cut_gpu_memory.h ljcl_cut_gpu.cpp
	$(OCL) -o $@ -c ljcl_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_cut_gpu_cl.h: lj96_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh lj96_cut_gpu_kernel.cu $(OBJ_DIR)/lj96_cut_gpu_cl.h;

$(OBJ_DIR)/lj96_cut_gpu_memory.o: $(ALL_H) lj96_cut_gpu_memory.h lj96_cut_gpu_memory.cpp  $(OBJ_DIR)/lj96_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/lj96_cut_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c lj96_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_cut_gpu.o: $(ALL_H) lj96_cut_gpu_memory.h lj96_cut_gpu.cpp
	$(OCL) -o $@ -c lj96_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmm_cut_gpu_cl.h: cmm_cut_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh cmm_cut_gpu_kernel.cu $(OBJ_DIR)/cmm_cut_gpu_cl.h;

$(OBJ_DIR)/cmm_cut_gpu_memory.o: $(ALL_H) cmm_cut_gpu_memory.h cmm_cut_gpu_memory.cpp  $(OBJ_DIR)/cmm_cut_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/cmm_cut_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c cmm_cut_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmm_cut_gpu.o: $(ALL_H) cmm_cut_gpu_memory.h cmm_cut_gpu.cpp
	$(OCL) -o $@ -c cmm_cut_gpu.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmmc_long_gpu_cl.h: cmmc_long_gpu_kernel.cu
	$(BSH) ./geryon/file_to_cstr.sh cmmc_long_gpu_kernel.cu $(OBJ_DIR)/cmmc_long_gpu_cl.h;

$(OBJ_DIR)/cmmc_long_gpu_memory.o: $(ALL_H) cmmc_long_gpu_memory.h cmmc_long_gpu_memory.cpp  $(OBJ_DIR)/cmmc_long_gpu_cl.h $(OBJ_DIR)/pair_gpu_nbor_cl.h $(OBJ_DIR)/cmmc_long_gpu_cl.h $(OBJ_DIR)/atomic_gpu_memory.o
	$(OCL) -o $@ -c cmmc_long_gpu_memory.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cmmc_long_gpu.o: $(ALL_H) cmmc_long_gpu_memory.h cmmc_long_gpu.cpp
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

