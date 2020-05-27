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
       $(OBJ_DIR)/lal_base_dipole.o $(OBJ_DIR)/lal_base_three.o \
       $(OBJ_DIR)/lal_base_dpd.o \
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
       $(OBJ_DIR)/lal_lj_sdk.o $(OBJ_DIR)/lal_lj_sdk_ext.o \
       $(OBJ_DIR)/lal_lj_sdk_long.o $(OBJ_DIR)/lal_lj_sdk_long_ext.o \
       $(OBJ_DIR)/lal_eam.o $(OBJ_DIR)/lal_eam_ext.o \
       $(OBJ_DIR)/lal_eam_fs_ext.o $(OBJ_DIR)/lal_eam_alloy_ext.o \
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
       $(OBJ_DIR)/lal_coul_dsf.o $(OBJ_DIR)/lal_coul_dsf_ext.o \
       $(OBJ_DIR)/lal_sw.o $(OBJ_DIR)/lal_sw_ext.o \
       $(OBJ_DIR)/lal_vashishta.o $(OBJ_DIR)/lal_vashishta_ext.o \
       $(OBJ_DIR)/lal_beck.o $(OBJ_DIR)/lal_beck_ext.o \
       $(OBJ_DIR)/lal_mie.o $(OBJ_DIR)/lal_mie_ext.o \
       $(OBJ_DIR)/lal_soft.o $(OBJ_DIR)/lal_soft_ext.o \
       $(OBJ_DIR)/lal_lj_coul_msm.o $(OBJ_DIR)/lal_lj_coul_msm_ext.o \
       $(OBJ_DIR)/lal_lj_gromacs.o $(OBJ_DIR)/lal_lj_gromacs_ext.o \
       $(OBJ_DIR)/lal_dpd.o $(OBJ_DIR)/lal_dpd_ext.o \
       $(OBJ_DIR)/lal_tersoff.o $(OBJ_DIR)/lal_tersoff_ext.o \
       $(OBJ_DIR)/lal_tersoff_zbl.o $(OBJ_DIR)/lal_tersoff_zbl_ext.o \
       $(OBJ_DIR)/lal_tersoff_mod.o $(OBJ_DIR)/lal_tersoff_mod_ext.o \
       $(OBJ_DIR)/lal_coul.o $(OBJ_DIR)/lal_coul_ext.o \
       $(OBJ_DIR)/lal_coul_debye.o $(OBJ_DIR)/lal_coul_debye_ext.o \
       $(OBJ_DIR)/lal_zbl.o $(OBJ_DIR)/lal_zbl_ext.o \
       $(OBJ_DIR)/lal_lj_cubic.o $(OBJ_DIR)/lal_lj_cubic_ext.o \
       $(OBJ_DIR)/lal_ufm.o $(OBJ_DIR)/lal_ufm_ext.o \
       $(OBJ_DIR)/lal_dipole_long_lj.o $(OBJ_DIR)/lal_dipole_long_lj_ext.o \
       $(OBJ_DIR)/lal_lj_expand_coul_long.o $(OBJ_DIR)/lal_lj_expand_coul_long_ext.o \
       $(OBJ_DIR)/lal_coul_long_cs.o $(OBJ_DIR)/lal_coul_long_cs_ext.o \
       $(OBJ_DIR)/lal_born_coul_long_cs.o $(OBJ_DIR)/lal_born_coul_long_cs_ext.o \
       $(OBJ_DIR)/lal_born_coul_wolf_cs.o $(OBJ_DIR)/lal_born_coul_wolf_cs_ext.o \
       $(OBJ_DIR)/lal_lj_tip4p_long.o $(OBJ_DIR)/lal_lj_tip4p_long_ext.o

KERS = $(OBJ_DIR)/device_cl.h $(OBJ_DIR)/atom_cl.h \
       $(OBJ_DIR)/neighbor_cpu_cl.h $(OBJ_DIR)/pppm_cl.h \
       $(OBJ_DIR)/ellipsoid_nbor_cl.h $(OBJ_DIR)/gayberne_cl.h \
       $(OBJ_DIR)/gayberne_lj_cl.h $(OBJ_DIR)/re_squared_cl.h \
       $(OBJ_DIR)/re_squared_lj_cl.h $(OBJ_DIR)/lj_cl.h $(OBJ_DIR)/lj96_cl.h \
       $(OBJ_DIR)/lj_expand_cl.h $(OBJ_DIR)/lj_coul_cl.h \
       $(OBJ_DIR)/lj_coul_long_cl.h $(OBJ_DIR)/lj_dsf_cl.h \
       $(OBJ_DIR)/lj_class2_long_cl.h \
       $(OBJ_DIR)/coul_long_cl.h $(OBJ_DIR)/morse_cl.h \
       $(OBJ_DIR)/charmm_long_cl.h $(OBJ_DIR)/lj_sdk_cl.h \
       $(OBJ_DIR)/lj_sdk_long_cl.h $(OBJ_DIR)/neighbor_gpu_cl.h \
       $(OBJ_DIR)/eam_cl.h $(OBJ_DIR)/buck_cl.h \
       $(OBJ_DIR)/buck_coul_cl.h $(OBJ_DIR)/buck_coul_long_cl.h \
       $(OBJ_DIR)/table_cl.h $(OBJ_DIR)/yukawa_cl.h \
       $(OBJ_DIR)/born_cl.h $(OBJ_DIR)/born_coul_wolf_cl.h \
       $(OBJ_DIR)/born_coul_long_cl.h $(OBJ_DIR)/dipole_lj_cl.h \
       $(OBJ_DIR)/dipole_lj_sf_cl.h $(OBJ_DIR)/colloid_cl.h \
       $(OBJ_DIR)/gauss_cl.h $(OBJ_DIR)/yukawa_colloid_cl.h \
       $(OBJ_DIR)/lj_coul_debye_cl.h $(OBJ_DIR)/coul_dsf_cl.h \
       $(OBJ_DIR)/sw_cl.h $(OBJ_DIR)/beck_cl.h $(OBJ_DIR)/mie_cl.h \
       $(OBJ_DIR)/soft_cl.h $(OBJ_DIR)/lj_coul_msm_cl.h \
       $(OBJ_DIR)/lj_gromacs_cl.h $(OBJ_DIR)/dpd_cl.h \
       $(OBJ_DIR)/lj_gauss_cl.h $(OBJ_DIR)/dzugutov_cl.h \
       $(OBJ_DIR)/tersoff_cl.h $(OBJ_DIR)/tersoff_zbl_cl.h \
       $(OBJ_DIR)/tersoff_mod_cl.h $(OBJ_DIR)/coul_cl.h \
       $(OBJ_DIR)/coul_debye_cl.h $(OBJ_DIR)/zbl_cl.h \
       $(OBJ_DIR)/lj_cubic_cl.h $(OBJ_DIR)/vashishta_cl.h \
       $(OBJ_DIR)/ufm_cl.h  $(OBJ_DIR)/dipole_long_lj_cl.h \
       $(OBJ_DIR)/lj_expand_coul_long_cl.h $(OBJ_DIR)/coul_long_cs_cl.h \
       $(OBJ_DIR)/born_coul_long_cs_cl.h $(OBJ_DIR)/born_coul_wolf_cs_cl.h \
       $(OBJ_DIR)/lj_tip4p_long_cl.h


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

$(OBJ_DIR)/lal_base_three.o: $(OCL_H) lal_base_three.h lal_base_three.cpp
	$(OCL) -o $@ -c lal_base_three.cpp

$(OBJ_DIR)/lal_base_dpd.o: $(OCL_H) lal_base_dpd.h lal_base_dpd.cpp
	$(OCL) -o $@ -c lal_base_dpd.cpp

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

$(OBJ_DIR)/lj_tip4p_long_cl.h: lal_lj_tip4p_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_tip4p_long $(PRE1_H) lal_lj_tip4p_long.cu $(OBJ_DIR)/lj_tip4p_long_cl.h;

$(OBJ_DIR)/lal_lj_tip4p_long.o: $(ALL_H) lal_lj_tip4p_long.h lal_lj_tip4p_long.cpp  $(OBJ_DIR)/lj_tip4p_long_cl.h $(OBJ_DIR)/lj_tip4p_long_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_lj_tip4p_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_tip4p_long_ext.o: $(ALL_H) lal_lj_tip4p_long.h lal_lj_tip4p_long_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_lj_tip4p_long_ext.cpp -I$(OBJ_DIR)

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

$(OBJ_DIR)/lj_sdk_cl.h: lal_lj_sdk.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_sdk $(PRE1_H) lal_lj_sdk.cu $(OBJ_DIR)/lj_sdk_cl.h;

$(OBJ_DIR)/lal_lj_sdk.o: $(ALL_H) lal_lj_sdk.h lal_lj_sdk.cpp  $(OBJ_DIR)/lj_sdk_cl.h $(OBJ_DIR)/lj_sdk_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_lj_sdk.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_sdk_ext.o: $(ALL_H) lal_lj_sdk.h lal_lj_sdk_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_lj_sdk_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_sdk_long_cl.h: lal_lj_sdk_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_sdk_long $(PRE1_H) lal_lj_sdk_long.cu $(OBJ_DIR)/lj_sdk_long_cl.h;

$(OBJ_DIR)/lal_lj_sdk_long.o: $(ALL_H) lal_lj_sdk_long.h lal_lj_sdk_long.cpp  $(OBJ_DIR)/lj_sdk_long_cl.h $(OBJ_DIR)/lj_sdk_long_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_lj_sdk_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_sdk_long_ext.o: $(ALL_H) lal_lj_sdk_long.h lal_lj_sdk_long_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_lj_sdk_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/eam_cl.h: lal_eam.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh eam $(PRE1_H) lal_eam.cu $(OBJ_DIR)/eam_cl.h;

$(OBJ_DIR)/lal_eam.o: $(ALL_H) lal_eam.h lal_eam.cpp  $(OBJ_DIR)/eam_cl.h $(OBJ_DIR)/eam_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_eam.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_ext.o: $(ALL_H) lal_eam.h lal_eam_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_eam_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_fs_ext.o: $(ALL_H) lal_eam.h lal_eam_fs_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_eam_fs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_alloy_ext.o: $(ALL_H) lal_eam.h lal_eam_alloy_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_eam_alloy_ext.cpp -I$(OBJ_DIR)

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

$(OBJ_DIR)/sw_cl.h: lal_sw.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh sw $(PRE1_H) lal_sw.cu $(OBJ_DIR)/sw_cl.h;

$(OBJ_DIR)/lal_sw.o: $(ALL_H) lal_sw.h lal_sw.cpp  $(OBJ_DIR)/sw_cl.h $(OBJ_DIR)/sw_cl.h $(OBJ_DIR)/lal_base_three.o
	$(OCL) -o $@ -c lal_sw.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_sw_ext.o: $(ALL_H) lal_sw.h lal_sw_ext.cpp lal_base_three.h
	$(OCL) -o $@ -c lal_sw_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/vashishta_cl.h: lal_vashishta.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh vashishta $(PRE1_H) lal_vashishta.cu $(OBJ_DIR)/vashishta_cl.h;

$(OBJ_DIR)/lal_vashishta.o: $(ALL_H) lal_vashishta.h lal_vashishta.cpp  $(OBJ_DIR)/vashishta_cl.h $(OBJ_DIR)/vashishta_cl.h $(OBJ_DIR)/lal_base_three.o
	$(OCL) -o $@ -c lal_vashishta.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_vashishta_ext.o: $(ALL_H) lal_vashishta.h lal_vashishta_ext.cpp lal_base_three.h
	$(OCL) -o $@ -c lal_vashishta_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/beck_cl.h: lal_beck.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh beck $(PRE1_H) lal_beck.cu $(OBJ_DIR)/beck_cl.h;

$(OBJ_DIR)/lal_beck.o: $(ALL_H) lal_beck.h lal_beck.cpp  $(OBJ_DIR)/beck_cl.h $(OBJ_DIR)/beck_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_beck.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_beck_ext.o: $(ALL_H) lal_beck.h lal_beck_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_beck_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/mie_cl.h: lal_mie.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh mie $(PRE1_H) lal_mie.cu $(OBJ_DIR)/mie_cl.h;

$(OBJ_DIR)/lal_mie.o: $(ALL_H) lal_mie.h lal_mie.cpp  $(OBJ_DIR)/mie_cl.h $(OBJ_DIR)/mie_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_mie.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_mie_ext.o: $(ALL_H) lal_mie.h lal_mie_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_mie_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/soft_cl.h: lal_soft.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh soft $(PRE1_H) lal_soft.cu $(OBJ_DIR)/soft_cl.h;

$(OBJ_DIR)/lal_soft.o: $(ALL_H) lal_soft.h lal_soft.cpp  $(OBJ_DIR)/soft_cl.h $(OBJ_DIR)/soft_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_soft.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_soft_ext.o: $(ALL_H) lal_soft.h lal_soft_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_soft_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_msm_cl.h: lal_lj_coul_msm.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_coul_msm $(PRE1_H) lal_lj_coul_msm.cu $(OBJ_DIR)/lj_coul_msm_cl.h;

$(OBJ_DIR)/lal_lj_coul_msm.o: $(ALL_H) lal_lj_coul_msm.h lal_lj_coul_msm.cpp  $(OBJ_DIR)/lj_coul_msm_cl.h $(OBJ_DIR)/lj_coul_msm_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_lj_coul_msm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_msm_ext.o: $(ALL_H) lal_lj_coul_msm.h lal_lj_coul_msm_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_lj_coul_msm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_gromacs_cl.h: lal_lj_gromacs.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_gromacs $(PRE1_H) lal_lj_gromacs.cu $(OBJ_DIR)/lj_gromacs_cl.h;

$(OBJ_DIR)/lal_lj_gromacs.o: $(ALL_H) lal_lj_gromacs.h lal_lj_gromacs.cpp  $(OBJ_DIR)/lj_gromacs_cl.h $(OBJ_DIR)/lj_gromacs_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_lj_gromacs.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_gromacs_ext.o: $(ALL_H) lal_lj_gromacs.h lal_lj_gromacs_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_lj_gromacs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/dpd_cl.h: lal_dpd.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh dpd $(PRE1_H) lal_dpd.cu $(OBJ_DIR)/dpd_cl.h;

$(OBJ_DIR)/lal_dpd.o: $(ALL_H) lal_dpd.h lal_dpd.cpp  $(OBJ_DIR)/dpd_cl.h $(OBJ_DIR)/dpd_cl.h $(OBJ_DIR)/lal_base_dpd.o
	$(OCL) -o $@ -c lal_dpd.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dpd_ext.o: $(ALL_H) lal_dpd.h lal_dpd_ext.cpp lal_base_dpd.h
	$(OCL) -o $@ -c lal_dpd_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/tersoff_cl.h: lal_tersoff.cu lal_tersoff_extra.h $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh tersoff $(PRE1_H) lal_tersoff_extra.h lal_tersoff.cu $(OBJ_DIR)/tersoff_cl.h;

$(OBJ_DIR)/lal_tersoff.o: $(ALL_H) lal_tersoff.h lal_tersoff.cpp  $(OBJ_DIR)/tersoff_cl.h $(OBJ_DIR)/tersoff_cl.h $(OBJ_DIR)/lal_base_three.o
	$(OCL) -o $@ -c lal_tersoff.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_tersoff_ext.o: $(ALL_H) lal_tersoff.h lal_tersoff_ext.cpp lal_base_three.h
	$(OCL) -o $@ -c lal_tersoff_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/tersoff_zbl_cl.h: lal_tersoff_zbl.cu lal_tersoff_zbl_extra.h $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh tersoff_zbl $(PRE1_H) lal_tersoff_zbl_extra.h lal_tersoff_zbl.cu $(OBJ_DIR)/tersoff_zbl_cl.h;

$(OBJ_DIR)/lal_tersoff_zbl.o: $(ALL_H) lal_tersoff_zbl.h lal_tersoff_zbl.cpp  $(OBJ_DIR)/tersoff_zbl_cl.h $(OBJ_DIR)/tersoff_zbl_cl.h $(OBJ_DIR)/lal_base_three.o
	$(OCL) -o $@ -c lal_tersoff_zbl.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_tersoff_zbl_ext.o: $(ALL_H) lal_tersoff_zbl.h lal_tersoff_zbl_ext.cpp lal_base_three.h
	$(OCL) -o $@ -c lal_tersoff_zbl_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/tersoff_mod_cl.h: lal_tersoff_mod.cu lal_tersoff_mod_extra.h $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh tersoff_mod $(PRE1_H) lal_tersoff_mod_extra.h lal_tersoff_mod.cu $(OBJ_DIR)/tersoff_mod_cl.h;

$(OBJ_DIR)/lal_tersoff_mod.o: $(ALL_H) lal_tersoff_mod.h lal_tersoff_mod.cpp  $(OBJ_DIR)/tersoff_mod_cl.h $(OBJ_DIR)/tersoff_mod_cl.h $(OBJ_DIR)/lal_base_three.o
	$(OCL) -o $@ -c lal_tersoff_mod.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_tersoff_mod_ext.o: $(ALL_H) lal_tersoff_mod.h lal_tersoff_mod_ext.cpp lal_base_three.h
	$(OCL) -o $@ -c lal_tersoff_mod_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_cl.h: lal_coul.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh coul $(PRE1_H) lal_coul.cu $(OBJ_DIR)/coul_cl.h;

$(OBJ_DIR)/lal_coul.o: $(ALL_H) lal_coul.h lal_coul.cpp  $(OBJ_DIR)/coul_cl.h $(OBJ_DIR)/coul_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_ext.o: $(ALL_H) lal_coul.h lal_coul_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_debye_cl.h: lal_coul_debye.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh coul_debye $(PRE1_H) lal_coul_debye.cu $(OBJ_DIR)/coul_debye_cl.h;

$(OBJ_DIR)/lal_coul_debye.o: $(ALL_H) lal_coul_debye.h lal_coul_debye.cpp  $(OBJ_DIR)/coul_debye_cl.h $(OBJ_DIR)/coul_debye_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_coul_debye.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_debye_ext.o: $(ALL_H) lal_coul_debye.h lal_coul_debye_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_coul_debye_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/zbl_cl.h: lal_zbl.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh zbl $(PRE1_H) lal_zbl.cu $(OBJ_DIR)/zbl_cl.h;

$(OBJ_DIR)/lal_zbl.o: $(ALL_H) lal_zbl.h lal_zbl.cpp  $(OBJ_DIR)/zbl_cl.h $(OBJ_DIR)/zbl_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_zbl.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_zbl_ext.o: $(ALL_H) lal_zbl.h lal_zbl_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_zbl_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cubic_cl.h: lal_lj_cubic.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_cubic $(PRE1_H) lal_lj_cubic.cu $(OBJ_DIR)/lj_cubic_cl.h;

$(OBJ_DIR)/lal_lj_cubic.o: $(ALL_H) lal_lj_cubic.h lal_lj_cubic.cpp  $(OBJ_DIR)/lj_cubic_cl.h $(OBJ_DIR)/lj_cubic_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_lj_cubic.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_cubic_ext.o: $(ALL_H) lal_lj_cubic.h lal_lj_cubic_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_lj_cubic_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ufm_cl.h: lal_ufm.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh ufm $(PRE1_H) lal_ufm.cu $(OBJ_DIR)/ufm_cl.h;

$(OBJ_DIR)/lal_ufm.o: $(ALL_H) lal_ufm.h lal_ufm.cpp  $(OBJ_DIR)/ufm_cl.h $(OBJ_DIR)/ufm_cl.h $(OBJ_DIR)/lal_base_atomic.o
	$(OCL) -o $@ -c lal_ufm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_ufm_ext.o: $(ALL_H) lal_ufm.h lal_ufm_ext.cpp lal_base_atomic.h
	$(OCL) -o $@ -c lal_ufm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/dipole_long_lj_cl.h: lal_dipole_long_lj.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh dipole_long_lj $(PRE1_H) lal_dipole_long_lj.cu $(OBJ_DIR)/dipole_long_lj_cl.h;

$(OBJ_DIR)/lal_dipole_long_lj.o: $(ALL_H) lal_dipole_long_lj.h lal_dipole_long_lj.cpp  $(OBJ_DIR)/dipole_long_lj_cl.h $(OBJ_DIR)/lj_expand_coul_long_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_dipole_long_lj.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dipole_long_lj_ext.o: $(ALL_H) lal_dipole_long_lj.h lal_dipole_long_lj_ext.cpp lal_base_dipole.h
	$(OCL) -o $@ -c lal_dipole_long_lj_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_coul_long_cl.h: lal_lj_expand_coul_long.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh lj_expand_coul_long $(PRE1_H) lal_lj_expand_coul_long.cu $(OBJ_DIR)/lj_expand_coul_long_cl.h;

$(OBJ_DIR)/lal_lj_expand_coul_long.o: $(ALL_H) lal_lj_expand_coul_long.h lal_lj_expand_coul_long.cpp  $(OBJ_DIR)/lj_expand_coul_long_cl.h $(OBJ_DIR)/lj_expand_coul_long_cl.h $(OBJ_DIR)/lal_base_charge.o
	$(OCL) -o $@ -c lal_lj_expand_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_expand_coul_long_ext.o: $(ALL_H) lal_lj_expand_coul_long.h lal_lj_expand_coul_long_ext.cpp lal_base_charge.h
	$(OCL) -o $@ -c lal_lj_expand_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_long_cs_cl.h: lal_coul_long_cs.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh coul_long_cs $(PRE1_H) lal_coul_long_cs.cu $(OBJ_DIR)/coul_long_cs_cl.h;

$(OBJ_DIR)/lal_coul_long_cs.o: $(ALL_H) lal_coul_long_cs.h lal_coul_long_cs.cpp  $(OBJ_DIR)/coul_long_cs_cl.h $(OBJ_DIR)/coul_long_cs_cl.h $(OBJ_DIR)/lal_base_charge.o $(OBJ_DIR)/lal_coul_long.o
	$(OCL) -o $@ -c lal_coul_long_cs.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_long_cs_ext.o: $(ALL_H) lal_coul_long_cs.h lal_coul_long_cs_ext.cpp lal_coul_long.h
	$(OCL) -o $@ -c lal_coul_long_cs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born_coul_long_cs_cl.h: lal_born_coul_long_cs.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh born_coul_long_cs $(PRE1_H) lal_born_coul_long_cs.cu $(OBJ_DIR)/born_coul_long_cs_cl.h;

$(OBJ_DIR)/lal_born_coul_long_cs.o: $(ALL_H) lal_born_coul_long_cs.h lal_born_coul_long_cs.cpp  $(OBJ_DIR)/born_coul_long_cs_cl.h $(OBJ_DIR)/born_coul_long_cs_cl.h $(OBJ_DIR)/lal_base_charge.o $(OBJ_DIR)/lal_born_coul_long.o
	$(OCL) -o $@ -c lal_born_coul_long_cs.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_coul_long_cs_ext.o: $(ALL_H) lal_born_coul_long_cs.h lal_born_coul_long_cs_ext.cpp lal_born_coul_long.h
	$(OCL) -o $@ -c lal_born_coul_long_cs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born_coul_wolf_cs_cl.h: lal_born_coul_wolf_cs.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh born_coul_wolf_cs $(PRE1_H) lal_born_coul_wolf_cs.cu $(OBJ_DIR)/born_coul_wolf_cs_cl.h;

$(OBJ_DIR)/lal_born_coul_wolf_cs.o: $(ALL_H) lal_born_coul_wolf_cs.h lal_born_coul_wolf_cs.cpp  $(OBJ_DIR)/born_coul_wolf_cs_cl.h $(OBJ_DIR)/born_coul_wolf_cs_cl.h $(OBJ_DIR)/lal_base_charge.o $(OBJ_DIR)/lal_born_coul_wolf.o
	$(OCL) -o $@ -c lal_born_coul_wolf_cs.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_coul_wolf_cs_ext.o: $(ALL_H) lal_born_coul_wolf_cs.h lal_born_coul_wolf_cs_ext.cpp lal_born_coul_wolf.h
	$(OCL) -o $@ -c lal_born_coul_wolf_cs_ext.cpp -I$(OBJ_DIR)

$(BIN_DIR)/ocl_get_devices: ./geryon/ucl_get_devices.cpp $(OCL_H)
	$(OCL) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_OPENCL $(OCL_LINK) 

$(OCL_LIB): $(OBJS) $(PTXS)
	$(AR) -crusv $(OCL_LIB) $(OBJS)
	@cp $(EXTRAMAKE) Makefile.lammps

opencl: $(OCL_EXECS)

clean:
	-rm -rf $(EXECS) $(OCL_EXECS) $(OCL_LIB) $(OBJS) $(KERS) *.linkinfo

veryclean: clean
	-rm -rf *~ *.linkinfo

