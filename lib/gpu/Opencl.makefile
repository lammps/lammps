OCL  = $(OCL_CPP) $(OCL_PREC) $(OCL_TUNE) -DUSE_OPENCL
OCL_LIB = $(LIB_DIR)/libgpu.a
# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
OCL_H  = $(wildcard ./geryon/ocl*.h) $(UCL_H)
# Headers for Pair Stuff
PAIR_H  = lal_atom.h lal_answer.h lal_neighbor_shared.h \
          lal_neighbor.h lal_precision.h lal_device.h \
          lal_balance.h lal_pppm.h
# Headers for Preprocessor/Auxiliary Functions
PRE1_H = lal_preprocessor.h lal_aux_fun1.h

ALL_H = $(OCL_H) $(PAIR_H)

EXECS = $(BIN_DIR)/ocl_get_devices
OBJS = $(OBJ_DIR)/lal_atom.o $(OBJ_DIR)/lal_answer.o \
       $(OBJ_DIR)/lal_neighbor_shared.o $(OBJ_DIR)/lal_neighbor.o \
       $(OBJ_DIR)/lal_device.o $(OBJ_DIR)/lal_base_atomic.o \
       $(OBJ_DIR)/lal_base_charge.o $(OBJ_DIR)/lal_base_ellipsoid.o \
       $(OBJ_DIR)/lal_base_dipole.o \
       $(OBJ_DIR)/lal_pppm.o $(OBJ_DIR)/lal_pppm_ext.o \
       $(OBJ_DIR)/lal_gayberne.o $(OBJ_DIR)/lal_gayberne_ext.o \
       $(OBJ_DIR)/lal_re_squared.o $(OBJ_DIR)/lal_re_squared_ext.o \
       $(OBJ_DIR)/lal_lj.o $(OBJ_DIR)/lal_lj_ext.o \
       $(OBJ_DIR)/lal_lj96.o $(OBJ_DIR)/lal_lj96_ext.o \
       $(OBJ_DIR)/lal_lj_expand.o $(OBJ_DIR)/lal_lj_expand_ext.o \
       $(OBJ_DIR)/lal_lj_coul.o $(OBJ_DIR)/lal_lj_coul_ext.o \
       $(OBJ_DIR)/lal_lj_coul_long.o $(OBJ_DIR)/lal_lj_coul_long_ext.o \
       $(OBJ_DIR)/lal_lj_dsf.o $(OBJ_DIR)/lal_lj_dsf_ext.o \
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
       $(OBJ_DIR)/lal_yukawa.o $(OBJ_DIR)/lal_yukawa_ext.o \
       $(OBJ_DIR)/lal_born.o $(OBJ_DIR)/lal_born_ext.o \
       $(OBJ_DIR)/lal_born_coul_wolf.o $(OBJ_DIR)/lal_born_coul_wolf_ext.o \
       $(OBJ_DIR)/lal_born_coul_long.o $(OBJ_DIR)/lal_born_coul_long_ext.o \
       $(OBJ_DIR)/lal_dipole_lj.o $(OBJ_DIR)/lal_dipole_lj_ext.o \
       $(OBJ_DIR)/lal_dipole_lj_sf.o $(OBJ_DIR)/lal_dipole_lj_sf_ext.o \
       $(OBJ_DIR)/lal_colloid.o $(OBJ_DIR)/lal_colloid_ext.o \
       $(OBJ_DIR)/lal_gauss.o $(OBJ_DIR)/lal_gauss_ext.o \
       $(OBJ_DIR)/lal_yukawa_colloid.o $(OBJ_DIR)/lal_yukawa_colloid_ext.o \
       $(OBJ_DIR)/lal_lj_coul_debye.o $(OBJ_DIR)/lal_lj_coul_debye_ext.o \
       $(OBJ_DIR)/lal_coul_dsf.o $(OBJ_DIR)/lal_coul_dsf_ext.o

KERS = $(OBJ_DIR)/device_cl.h $(OBJ_DIR)/atom_cl.h \
       $(OBJ_DIR)/neighbor_cpu_cl.h $(OBJ_DIR)/pppm_cl.h \
       $(OBJ_DIR)/ellipsoid_nbor_cl.h $(OBJ_DIR)/gayberne_cl.h \
       $(OBJ_DIR)/gayberne_lj_cl.h $(OBJ_DIR)/re_squared_cl.h \
       $(OBJ_DIR)/re_squared_lj_cl.h $(OBJ_DIR)/lj_cl.h $(OBJ_DIR)/lj96_cl.h \
       $(OBJ_DIR)/lj_expand_cl.h $(OBJ_DIR)/lj_coul_cl.h \
       $(OBJ_DIR)/lj_coul_long_cl.h $(OBJ_DIR)/lj_dsf_cl.h \
       $(OBJ_DIR)/lj_class2_long_cl.h \
       $(OBJ_DIR)/coul_long_cl.h $(OBJ_DIR)/morse_cl.h \
       $(OBJ_DIR)/charmm_long_cl.h $(OBJ_DIR)/cg_cmm_cl.h \
       $(OBJ_DIR)/cg_cmm_long_cl.h $(OBJ_DIR)/neighbor_gpu_cl.h \
       $(OBJ_DIR)/eam_cl.h $(OBJ_DIR)/buck_cl.h \
       $(OBJ_DIR)/buck_coul_cl.h $(OBJ_DIR)/buck_coul_long_cl.h \
       $(OBJ_DIR)/table_cl.h $(OBJ_DIR)/yukawa_cl.h \
       $(OBJ_DIR)/born_cl.h $(OBJ_DIR)/born_coul_wolf_cl.h \
       $(OBJ_DIR)/born_coul_long_cl.h $(OBJ_DIR)/dipole_lj_cl.h \
       $(OBJ_DIR)/dipole_lj_sf_cl.h $(OBJ_DIR)/colloid_cl.h \
       $(OBJ_DIR)/gauss_cl.h $(OBJ_DIR)/yukawa_colloid_cl.h \
       $(OBJ_DIR)/lj_coul_debye_cl.h $(OBJ_DIR)/coul_dsf_cl.h


OCL_EXECS = $(BIN_DIR)/ocl_get_devices

all: $(OBJ_DIR) $(OCL_LIB) $(EXECS)

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/atom_cl.h: lal_atom.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh atom lal_preprocessor.h lal_atom.cu $(OBJ_DIR)/atom_cl.h

$(OBJ_DIR)/lal_atom.o: lal_atom.cpp lal_atom.h $(OCL_H) $(OBJ_DIR)/atom_cl.h
	$(OCL) -o $@ -c lal_atom.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_answer.o: lal_answer.cpp lal_answer.h $(OCL_H)
	$(OCL) -o $@ -c lal_answer.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/neighbor_cpu_cl.h: lal_neighbor_cpu.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh neighbor_cpu lal_preprocessor.h lal_neighbor_cpu.cu $(OBJ_DIR)/neighbor_cpu_cl.h

$(OBJ_DIR)/neighbor_gpu_cl.h: lal_neighbor_gpu.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh neighbor_gpu lal_preprocessor.h lal_neighbor_gpu.cu $(OBJ_DIR)/neighbor_gpu_cl.h

$(OBJ_DIR)/lal_neighbor_shared.o: lal_neighbor_shared.cpp lal_neighbor_shared.h $(OCL_H) $(OBJ_DIR)/neighbor_cpu_cl.h $(OBJ_DIR)/neighbor_gpu_cl.h
	$(OCL) -o $@ -c lal_neighbor_shared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_neighbor.o: lal_neighbor.cpp lal_neighbor.h $(OCL_H) lal_neighbor_shared.h
	$(OCL) -o $@ -c lal_neighbor.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/device_cl.h: lal_device.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh device lal_preprocessor.h lal_device.cu $(OBJ_DIR)/device_cl.h

$(OBJ_DIR)/lal_device.o: lal_device.cpp lal_device.h $(ALL_H) $(OBJ_DIR)/device_cl.h
	$(OCL) -o $@ -c lal_device.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_base_atomic.o: $(OCL_H) lal_base_atomic.h lal_base_atomic.cpp
	$(OCL) -o $@ -c lal_base_atomic.cpp

$(OBJ_DIR)/lal_base_charge.o: $(OCL_H) lal_base_charge.h lal_base_charge.cpp
	$(OCL) -o $@ -c lal_base_charge.cpp

$(OBJ_DIR)/lal_base_ellipsoid.o: $(OCL_H) lal_base_ellipsoid.h lal_base_ellipsoid.cpp $(OBJ_DIR)/ellipsoid_nbor_cl.h
	$(OCL) -o $@ -c lal_base_ellipsoid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_base_dipole.o: $(OCL_H) lal_base_dipole.h lal_base_dipole.cpp
	$(OCL) -o $@ -c lal_base_dipole.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pppm_cl.h: lal_pppm.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh pppm lal_preprocessor.h lal_pppm.cu $(OBJ_DIR)/pppm_cl.h;

$(OBJ_DIR)/lal_pppm.o: $(ALL_H) lal_pppm.h lal_pppm.cpp  $(OBJ_DIR)/pppm_cl.h $(OBJ_DIR)/pppm_cl.h
	$(OCL) -o $@ -c lal_pppm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_pppm_ext.o: $(ALL_H) lal_pppm.h lal_pppm_ext.cpp
	$(OCL) -o $@ -c lal_pppm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ellipsoid_nbor_cl.h: lal_ellipsoid_nbor.cu lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh ellipsoid_nbor lal_preprocessor.h lal_ellipsoid_nbor.cu $(OBJ_DIR)/ellipsoid_nbor_cl.h

$(OBJ_DIR)/gayberne_cl.h: lal_gayberne.cu lal_ellipsoid_extra.h lal_aux_fun1.h lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh gayberne lal_preprocessor.h lal_aux_fun1.h lal_ellipsoid_extra.h lal_gayberne.cu $(OBJ_DIR)/gayberne_cl.h;

$(OBJ_DIR)/gayberne_lj_cl.h: lal_gayberne_lj.cu lal_ellipsoid_extra.h lal_aux_fun1.h lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh gayberne_lj lal_preprocessor.h lal_aux_fun1.h lal_ellipsoid_extra.h lal_gayberne_lj.cu $(OBJ_DIR)/gayberne_lj_cl.h;

$(OBJ_DIR)/lal_gayberne.o: $(ALL_H) lal_gayberne.h lal_gayberne.cpp $(OBJ_DIR)/gayberne_cl.h $(OBJ_DIR)/gayberne_lj_cl.h $(OBJ_DIR)/lal_base_ellipsoid.o
	$(OCL) -o $@ -c lal_gayberne.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_gayberne_ext.o: $(ALL_H) $(OBJ_DIR)/lal_gayberne.o lal_gayberne_ext.cpp
	$(OCL) -o $@ -c lal_gayberne_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/re_squared_cl.h: lal_re_squared.cu lal_ellipsoid_extra.h lal_aux_fun1.h lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh re_squared lal_preprocessor.h lal_aux_fun1.h lal_ellipsoid_extra.h lal_re_squared.cu $(OBJ_DIR)/re_squared_cl.h;

$(OBJ_DIR)/re_squared_lj_cl.h: lal_re_squared_lj.cu lal_ellipsoid_extra.h lal_aux_fun1.h lal_preprocessor.h
	$(BSH) ./geryon/file_to_cstr.sh re_squared_lj lal_preprocessor.h lal_aux_fun1.h lal_ellipsoid_extra.h lal_re_squared_lj.cu $(OBJ_DIR)/re_squared_lj_cl.h;

$(OBJ_DIR)/lal_re_squared.o: $(ALL_H) lal_re_squared.h lal_re_squared.cpp $(OBJ_DIR)/re_squared_cl.h $(OBJ_DIR)/re_squared_lj_cl.h $(OBJ_DIR)/lal_base_ellipsoid.o
	$(OCL) -o $@ -c lal_re_squared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_re_squared_ext.o: $(ALL_H) $(OBJ_DIR)/lal_re_squared.o lal_re_squared_ext.cpp
	$(OCL) -o $@ -c lal_re_squared_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cl.h: lal_lj.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj $(PRE1_H) lal_lj.cu $(OBJ_DIR)/lj_cl.h;

$(OBJ_DIR)/lal_lj.o: $(ALL_H) lal_lj.h lal_lj.cpp  $(OBJ_DIR)/lj_cl.h $(OBJ_DIR)/lj_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_lj.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_ext.o: $(ALL_H) lal_lj.h lal_lj_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_lj_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_cl.h: lal_lj_coul.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_coul $(PRE1_H) lal_lj_coul.cu $(OBJ_DIR)/lj_coul_cl.h;

$(OBJ_DIR)/lal_lj_coul.o: $(ALL_H) lal_lj_coul.h lal_lj_coul.cpp  $(OBJ_DIR)/lj_coul_cl.h $(OBJ_DIR)/lj_coul_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_lj_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_ext.o: $(ALL_H) lal_lj_coul.h lal_lj_coul_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_lj_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_long_cl.h: lal_lj_coul_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_coul_long $(PRE1_H) lal_lj_coul_long.cu $(OBJ_DIR)/lj_coul_long_cl.h;

$(OBJ_DIR)/lal_lj_coul_long.o: $(ALL_H) lal_lj_coul_long.h lal_lj_coul_long.cpp  $(OBJ_DIR)/lj_coul_long_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_lj_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_long_ext.o: $(ALL_H) lal_lj_coul_long.h lal_lj_coul_long_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_lj_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_dsf_cl.h: lal_lj_dsf.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_dsf $(PRE1_H) lal_lj_dsf.cu $(OBJ_DIR)/lj_dsf_cl.h;

$(OBJ_DIR)/lal_lj_dsf.o: $(ALL_H) lal_lj_dsf.h lal_lj_dsf.cpp  $(OBJ_DIR)/lj_dsf_cl.h $(OBJ_DIR)/lj_dsf_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_lj_dsf.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_dsf_ext.o: $(ALL_H) lal_lj_dsf.h lal_lj_dsf_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_lj_dsf_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long_cl.h: lal_lj_class2_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_class2_long $(PRE1_H) lal_lj_class2_long.cu $(OBJ_DIR)/lj_class2_long_cl.h;

$(OBJ_DIR)/lal_lj_class2_long.o: $(ALL_H) lal_lj_class2_long.h lal_lj_class2_long.cpp  $(OBJ_DIR)/lj_class2_long_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_lj_class2_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_class2_long_ext.o: $(ALL_H) lal_lj_class2_long.h lal_lj_class2_long_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_lj_class2_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_long_cl.h: lal_coul_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh coul_long $(PRE1_H) lal_coul_long.cu $(OBJ_DIR)/coul_long_cl.h;

$(OBJ_DIR)/lal_coul_long.o: $(ALL_H) lal_coul_long.h lal_coul_long.cpp  $(OBJ_DIR)/coul_long_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_long_ext.o: $(ALL_H) lal_coul_long.h lal_coul_long_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse_cl.h: lal_morse.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh morse $(PRE1_H) lal_morse.cu $(OBJ_DIR)/morse_cl.h;

$(OBJ_DIR)/lal_morse.o: $(ALL_H) lal_morse.h lal_morse.cpp  $(OBJ_DIR)/morse_cl.h $(OBJ_DIR)/morse_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_morse.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_morse_ext.o: $(ALL_H) lal_morse.h lal_morse_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_morse_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/charmm_long_cl.h: lal_charmm_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh charmm_long $(PRE1_H) lal_charmm_long.cu $(OBJ_DIR)/charmm_long_cl.h;

$(OBJ_DIR)/lal_charmm_long.o: $(ALL_H) lal_charmm_long.h lal_charmm_long.cpp  $(OBJ_DIR)/charmm_long_cl.h $(OBJ_DIR)/charmm_long_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_charmm_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_charmm_long_ext.o: $(ALL_H) lal_charmm_long.h lal_charmm_long_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_charmm_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_cl.h: lal_lj96.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj96 $(PRE1_H) lal_lj96.cu $(OBJ_DIR)/lj96_cl.h;

$(OBJ_DIR)/lal_lj96.o: $(ALL_H) lal_lj96.h lal_lj96.cpp  $(OBJ_DIR)/lj96_cl.h $(OBJ_DIR)/lj96_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_lj96.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj96_ext.o: $(ALL_H) lal_lj96.h lal_lj96_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_lj96_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_cl.h: lal_lj_expand.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_expand $(PRE1_H) lal_lj_expand.cu $(OBJ_DIR)/lj_expand_cl.h;

$(OBJ_DIR)/lal_lj_expand.o: $(ALL_H) lal_lj_expand.h lal_lj_expand.cpp  $(OBJ_DIR)/lj_expand_cl.h $(OBJ_DIR)/lj_expand_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_lj_expand.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_expand_ext.o: $(ALL_H) lal_lj_expand.h lal_lj_expand_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_lj_expand_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_cl.h: lal_cg_cmm.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh cg_cmm $(PRE1_H) lal_cg_cmm.cu $(OBJ_DIR)/cg_cmm_cl.h;

$(OBJ_DIR)/lal_cg_cmm.o: $(ALL_H) lal_cg_cmm.h lal_cg_cmm.cpp  $(OBJ_DIR)/cg_cmm_cl.h $(OBJ_DIR)/cg_cmm_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_cg_cmm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_cg_cmm_ext.o: $(ALL_H) lal_cg_cmm.h lal_cg_cmm_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_cg_cmm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_long_cl.h: lal_cg_cmm_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh cg_cmm_long $(PRE1_H) lal_cg_cmm_long.cu $(OBJ_DIR)/cg_cmm_long_cl.h;

$(OBJ_DIR)/lal_cg_cmm_long.o: $(ALL_H) lal_cg_cmm_long.h lal_cg_cmm_long.cpp  $(OBJ_DIR)/cg_cmm_long_cl.h $(OBJ_DIR)/cg_cmm_long_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_cg_cmm_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_cg_cmm_long_ext.o: $(ALL_H) lal_cg_cmm_long.h lal_cg_cmm_long_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_cg_cmm_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/eam_cl.h: lal_eam.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh eam $(PRE1_H) lal_eam.cu $(OBJ_DIR)/eam_cl.h;

$(OBJ_DIR)/lal_eam.o: $(ALL_H) lal_eam.h lal_eam.cpp  $(OBJ_DIR)/eam_cl.h $(OBJ_DIR)/eam_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_eam.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_ext.o: $(ALL_H) lal_eam.h lal_eam_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_eam_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/buck_cl.h: lal_buck.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh buck $(PRE1_H) lal_buck.cu $(OBJ_DIR)/buck_cl.h;

$(OBJ_DIR)/lal_buck.o: $(ALL_H) lal_buck.h lal_buck.cpp  $(OBJ_DIR)/buck_cl.h $(OBJ_DIR)/buck_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_buck.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_buck_ext.o: $(ALL_H) lal_buck.h lal_buck_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_buck_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/buck_coul_cl.h: lal_buck_coul.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh buck_coul $(PRE1_H) lal_buck_coul.cu $(OBJ_DIR)/buck_coul_cl.h;

$(OBJ_DIR)/lal_buck_coul.o: $(ALL_H) lal_buck_coul.h lal_buck_coul.cpp  $(OBJ_DIR)/buck_coul_cl.h $(OBJ_DIR)/buck_coul_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_buck_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_buck_coul_ext.o: $(ALL_H) lal_buck_coul.h lal_buck_coul_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_buck_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/buck_coul_long_cl.h: lal_buck_coul_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh buck_coul_long $(PRE1_H) lal_buck_coul_long.cu $(OBJ_DIR)/buck_coul_long_cl.h;

$(OBJ_DIR)/lal_buck_coul_long.o: $(ALL_H) lal_buck_coul_long.h lal_buck_coul_long.cpp  $(OBJ_DIR)/buck_coul_long_cl.h $(OBJ_DIR)/buck_coul_long_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_buck_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_buck_coul_long_ext.o: $(ALL_H) lal_buck_coul_long.h lal_buck_coul_long_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_buck_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/table_cl.h: lal_table.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh table $(PRE1_H) lal_table.cu $(OBJ_DIR)/table_cl.h;

$(OBJ_DIR)/lal_table.o: $(ALL_H) lal_table.h lal_table.cpp  $(OBJ_DIR)/table_cl.h $(OBJ_DIR)/table_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_table.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_table_ext.o: $(ALL_H) lal_table.h lal_table_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_table_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/yukawa_cl.h: lal_yukawa.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh yukawa $(PRE1_H) lal_yukawa.cu $(OBJ_DIR)/yukawa_cl.h;

$(OBJ_DIR)/lal_yukawa.o: $(ALL_H) lal_yukawa.h lal_yukawa.cpp  $(OBJ_DIR)/yukawa_cl.h $(OBJ_DIR)/yukawa_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_yukawa.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_yukawa_ext.o: $(ALL_H) lal_yukawa.h lal_yukawa_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_yukawa_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born_cl.h: lal_born.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh born $(PRE1_H) lal_born.cu $(OBJ_DIR)/born_cl.h;

$(OBJ_DIR)/lal_born.o: $(ALL_H) lal_born.h lal_born.cpp  $(OBJ_DIR)/born_cl.h $(OBJ_DIR)/born_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_born.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_ext.o: $(ALL_H) lal_born.h lal_born_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_born_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born_coul_wolf_cl.h: lal_born_coul_wolf.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh born_coul_wolf $(PRE1_H) lal_born_coul_wolf.cu $(OBJ_DIR)/born_coul_wolf_cl.h;

$(OBJ_DIR)/lal_born_coul_wolf.o: $(ALL_H) lal_born_coul_wolf.h lal_born_coul_wolf.cpp  $(OBJ_DIR)/born_coul_wolf_cl.h $(OBJ_DIR)/born_coul_wolf_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_born_coul_wolf.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_coul_wolf_ext.o: $(ALL_H) lal_born_coul_wolf.h lal_born_coul_wolf_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_born_coul_wolf_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born_coul_long_cl.h: lal_born_coul_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh born_coul_long $(PRE1_H) lal_born_coul_long.cu $(OBJ_DIR)/born_coul_long_cl.h;

$(OBJ_DIR)/lal_born_coul_long.o: $(ALL_H) lal_born_coul_long.h lal_born_coul_long.cpp  $(OBJ_DIR)/born_coul_long_cl.h $(OBJ_DIR)/born_coul_long_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_born_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_coul_long_ext.o: $(ALL_H) lal_born_coul_long.h lal_born_coul_long_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_born_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/dipole_lj_cl.h: lal_dipole_lj.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh dipole_lj $(PRE1_H) lal_dipole_lj.cu $(OBJ_DIR)/dipole_lj_cl.h;

$(OBJ_DIR)/lal_dipole_lj.o: $(ALL_H) lal_dipole_lj.h lal_dipole_lj.cpp  $(OBJ_DIR)/dipole_lj_cl.h $(OBJ_DIR)/dipole_lj_cl.h $(OBJ_DIR)/lal_base_dipole.o
	$(OCL) -o $@ -c lal_dipole_lj.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dipole_lj_ext.o: $(ALL_H) lal_dipole_lj.h lal_dipole_lj_ext.cpp lal_base_dipole.h
	$(OCL) -o $@ -c lal_dipole_lj_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/dipole_lj_sf_cl.h: lal_dipole_lj_sf.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh dipole_lj_sf $(PRE1_H) lal_dipole_lj_sf.cu $(OBJ_DIR)/dipole_lj_sf_cl.h;

$(OBJ_DIR)/lal_dipole_lj_sf.o: $(ALL_H) lal_dipole_lj_sf.h lal_dipole_lj_sf.cpp  $(OBJ_DIR)/dipole_lj_sf_cl.h $(OBJ_DIR)/dipole_lj_sf_cl.h $(OBJ_DIR)/lal_base_dipole.o
	$(OCL) -o $@ -c lal_dipole_lj_sf.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dipole_lj_sf_ext.o: $(ALL_H) lal_dipole_lj_sf.h lal_dipole_lj_sf_ext.cpp lal_base_dipole.h
	$(OCL) -o $@ -c lal_dipole_lj_sf_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/colloid_cl.h: lal_colloid.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh colloid $(PRE1_H) lal_colloid.cu $(OBJ_DIR)/colloid_cl.h;

$(OBJ_DIR)/lal_colloid.o: $(ALL_H) lal_colloid.h lal_colloid.cpp  $(OBJ_DIR)/colloid_cl.h $(OBJ_DIR)/colloid_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_colloid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_colloid_ext.o: $(ALL_H) lal_colloid.h lal_colloid_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_colloid_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/gauss_cl.h: lal_gauss.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh gauss $(PRE1_H) lal_gauss.cu $(OBJ_DIR)/gauss_cl.h;

$(OBJ_DIR)/lal_gauss.o: $(ALL_H) lal_gauss.h lal_gauss.cpp  $(OBJ_DIR)/gauss_cl.h $(OBJ_DIR)/gauss_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_gauss.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_gauss_ext.o: $(ALL_H) lal_gauss.h lal_gauss_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_gauss_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/yukawa_colloid_cl.h: lal_yukawa_colloid.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh yukawa_colloid $(PRE1_H) lal_yukawa_colloid.cu $(OBJ_DIR)/yukawa_colloid_cl.h;

$(OBJ_DIR)/lal_yukawa_colloid.o: $(ALL_H) lal_yukawa_colloid.h lal_yukawa_colloid.cpp  $(OBJ_DIR)/yukawa_colloid_cl.h $(OBJ_DIR)/yukawa_colloid_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_yukawa_colloid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_yukawa_colloid_ext.o: $(ALL_H) lal_yukawa_colloid.h lal_yukawa_colloid_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_yukawa_colloid_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_debye_cl.h: lal_lj_coul_debye.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_coul_debye $(PRE1_H) lal_lj_coul_debye.cu $(OBJ_DIR)/lj_coul_debye_cl.h;

$(OBJ_DIR)/lal_lj_coul_debye.o: $(ALL_H) lal_lj_coul_debye.h lal_lj_coul_debye.cpp  $(OBJ_DIR)/lj_coul_debye_cl.h $(OBJ_DIR)/lj_coul_debye_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_lj_coul_debye.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_debye_ext.o: $(ALL_H) lal_lj_coul_debye.h lal_lj_coul_debye_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_lj_coul_debye_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_dsf_cl.h: lal_coul_dsf.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh coul_dsf $(PRE1_H) lal_coul_dsf.cu $(OBJ_DIR)/coul_dsf_cl.h;

$(OBJ_DIR)/lal_coul_dsf.o: $(ALL_H) lal_coul_dsf.h lal_coul_dsf.cpp  $(OBJ_DIR)/coul_dsf_cl.h $(OBJ_DIR)/coul_dsf_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_coul_dsf.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_dsf_ext.o: $(ALL_H) lal_coul_dsf.h lal_coul_dsf_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_coul_dsf_ext.cpp -I$(OBJ_DIR)

$(BIN_DIR)/ocl_get_devices: ./geryon/ucl_get_devices.cpp
	$(OCL) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_OPENCL $(OCL_LINK) 

$(OCL_LIB): $(OBJS) $(PTXS)
	$(AR) -crusv $(OCL_LIB) $(OBJS)
	@cp $(EXTRAMAKE) Makefile.lammps

opencl: $(OCL_EXECS)

clean:
	rm -rf $(EXECS) $(OCL_EXECS) $(OCL_LIB) $(OBJS) $(KERS) *.linkinfo

veryclean: clean
	rm -rf *~ *.linkinfo

