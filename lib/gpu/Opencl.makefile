OCL  = $(OCL_CPP) $(OCL_PREC) -DUSE_OPENCL
OCL_LIB = $(LIB_DIR)/libgpu.a
# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
OCL_H  = $(wildcard ./geryon/ocl*.h) $(UCL_H)
# Headers for Pair Stuff
PAIR_H  = atom.h answer.h neighbor_shared.h \
          neighbor.h precision.h device.h \
          balance.h pppm.h

ALL_H = $(OCL_H) $(PAIR_H)

EXECS = $(BIN_DIR)/ocl_get_devices
OBJS = $(OBJ_DIR)/atom.o $(OBJ_DIR)/ans.o \
       $(OBJ_DIR)/neighbor_shared.o $(OBJ_DIR)/neighbor.o \
       $(OBJ_DIR)/device.o $(OBJ_DIR)/base_atomic.o \
       $(OBJ_DIR)/base_charge.o $(OBJ_DIR)/base_ellipsoid.o \
       $(OBJ_DIR)/pppm.o $(OBJ_DIR)/pppm_ext.o \
       $(OBJ_DIR)/gayberne.o $(OBJ_DIR)/gayberne_ext.o \
       $(OBJ_DIR)/re_squared.o $(OBJ_DIR)/re_squared_ext.o \
       $(OBJ_DIR)/lj.o $(OBJ_DIR)/lj_ext.o \
       $(OBJ_DIR)/lj96.o $(OBJ_DIR)/lj96_ext.o \
       $(OBJ_DIR)/lj_expand.o $(OBJ_DIR)/lj_expand_ext.o \
       $(OBJ_DIR)/lj_coul.o $(OBJ_DIR)/lj_coul_ext.o \
       $(OBJ_DIR)/lj_coul_long.o $(OBJ_DIR)/lj_coul_long_ext.o \
       $(OBJ_DIR)/lj_class2_long.o $(OBJ_DIR)/lj_class2_long_ext.o \
       $(OBJ_DIR)/coul_long.o $(OBJ_DIR)/coul_long_ext.o \
       $(OBJ_DIR)/morse.o $(OBJ_DIR)/morse_ext.o \
       $(OBJ_DIR)/charmm_long.o $(OBJ_DIR)/charmm_long_ext.o \
       $(OBJ_DIR)/cg_cmm.o $(OBJ_DIR)/cg_cmm_ext.o \
       $(OBJ_DIR)/cg_cmm_long.o $(OBJ_DIR)/cg_cmm_long_ext.o 
KERS = $(OBJ_DIR)/device_cl.h $(OBJ_DIR)/atom_cl.h \
       $(OBJ_DIR)/neighbor_cpu_cl.h $(OBJ_DIR)/pppm_cl.h \
       $(OBJ_DIR)/ellipsoid_nbor_cl.h $(OBJ_DIR)/gayberne_cl.h \
       $(OBJ_DIR)/gayberne_lj_cl.h $(OBJ_DIR)/re_squared_cl.h \
       $(OBJ_DIR)/re_squared_lj_cl.h $(OBJ_DIR)/lj_cl.h $(OBJ_DIR)/lj96_cl.h \
       $(OBJ_DIR)/lj_expand_cl.h $(OBJ_DIR)/lj_coul_cl.h \
       $(OBJ_DIR)/lj_coul_long_cl.h $(OBJ_DIR)/lj_class2_long_cl.h \
       $(OBJ_DIR)/coul_long_cl.h $(OBJ_DIR)/morse_cl.h \
       $(OBJ_DIR)/charmm_long_cl.h $(OBJ_DIR)/cg_cmm_cl.h \
       $(OBJ_DIR)/cg_cmm_long_cl.h 

OCL_EXECS = $(BIN_DIR)/ocl_get_devices

all: $(OCL_LIB) $(EXECS)

$(OBJ_DIR)/atom_cl.h: atom.cu
	$(BSH) ./geryon/file_to_cstr.sh atom atom.cu $(OBJ_DIR)/atom_cl.h

$(OBJ_DIR)/atom.o: atom.cpp atom.h $(OCL_H) $(OBJ_DIR)/atom_cl.h
	$(OCL) -o $@ -c atom.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ans.o: answer.cpp answer.h $(OCL_H)
	$(OCL) -o $@ -c answer.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/neighbor_cpu_cl.h: neighbor_cpu.cu
	$(BSH) ./geryon/file_to_cstr.sh neighbor_cpu neighbor_cpu.cu $(OBJ_DIR)/neighbor_cpu_cl.h

$(OBJ_DIR)/neighbor_shared.o: neighbor_shared.cpp neighbor_shared.h $(OCL_H) $(OBJ_DIR)/neighbor_cpu_cl.h
	$(OCL) -o $@ -c neighbor_shared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/neighbor.o: neighbor.cpp neighbor.h $(OCL_H) neighbor_shared.h
	$(OCL) -o $@ -c neighbor.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/device_cl.h: device.cu
	$(BSH) ./geryon/file_to_cstr.sh device device.cu $(OBJ_DIR)/device_cl.h

$(OBJ_DIR)/device.o: device.cpp device.h $(ALL_H) $(OBJ_DIR)/device_cl.h
	$(OCL) -o $@ -c device.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/base_atomic.o: $(OCL_H) base_atomic.h base_atomic.cpp
	$(OCL) -o $@ -c base_atomic.cpp

$(OBJ_DIR)/base_charge.o: $(OCL_H) base_charge.h base_charge.cpp
	$(OCL) -o $@ -c base_charge.cpp

$(OBJ_DIR)/base_ellipsoid.o: $(OCL_H) base_ellipsoid.h base_ellipsoid.cpp $(OBJ_DIR)/ellipsoid_nbor_cl.h
	$(OCL) -o $@ -c base_ellipsoid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pppm_cl.h: pppm.cu
	$(BSH) ./geryon/file_to_cstr.sh pppm pppm.cu $(OBJ_DIR)/pppm_cl.h;

$(OBJ_DIR)/pppm.o: $(ALL_H) pppm.h pppm.cpp  $(OBJ_DIR)/pppm_cl.h $(OBJ_DIR)/pppm_cl.h
	$(OCL) -o $@ -c pppm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/pppm_ext.o: $(ALL_H) pppm.h pppm_ext.cpp
	$(OCL) -o $@ -c pppm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ellipsoid_nbor_cl.h: ellipsoid_nbor.cu
	$(BSH) ./geryon/file_to_cstr.sh ellipsoid_nbor ellipsoid_nbor.cu $(OBJ_DIR)/ellipsoid_nbor_cl.h

$(OBJ_DIR)/gayberne_cl.h: gayberne.cu ellipsoid_extra.h
	$(BSH) ./geryon/file_to_cstr.sh gayberne ellipsoid_extra.h gayberne.cu $(OBJ_DIR)/gayberne_cl.h;

$(OBJ_DIR)/gayberne_lj_cl.h: gayberne_lj.cu ellipsoid_extra.h
	$(BSH) ./geryon/file_to_cstr.sh gayberne_lj ellipsoid_extra.h gayberne_lj.cu $(OBJ_DIR)/gayberne_lj_cl.h;

$(OBJ_DIR)/gayberne.o: $(ALL_H) gayberne.h gayberne.cpp $(OBJ_DIR)/gayberne_cl.h $(OBJ_DIR)/gayberne_lj_cl.h $(OBJ_DIR)/base_ellipsoid.o
	$(OCL) -o $@ -c gayberne.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/gayberne_ext.o: $(ALL_H) $(OBJ_DIR)/gayberne.o gayberne_ext.cpp
	$(OCL) -o $@ -c gayberne_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/re_squared_cl.h: re_squared.cu ellipsoid_extra.h
	$(BSH) ./geryon/file_to_cstr.sh re_squared ellipsoid_extra.h re_squared.cu $(OBJ_DIR)/re_squared_cl.h;

$(OBJ_DIR)/re_squared_lj_cl.h: re_squared_lj.cu ellipsoid_extra.h
	$(BSH) ./geryon/file_to_cstr.sh re_squared_lj ellipsoid_extra.h re_squared_lj.cu $(OBJ_DIR)/re_squared_lj_cl.h;

$(OBJ_DIR)/re_squared.o: $(ALL_H) re_squared.h re_squared.cpp $(OBJ_DIR)/re_squared_cl.h $(OBJ_DIR)/re_squared_lj_cl.h $(OBJ_DIR)/base_ellipsoid.o
	$(OCL) -o $@ -c re_squared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/re_squared_ext.o: $(ALL_H) $(OBJ_DIR)/re_squared.o re_squared_ext.cpp
	$(OCL) -o $@ -c re_squared_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cl.h: lj.cu
	$(BSH) ./geryon/file_to_cstr.sh lj lj.cu $(OBJ_DIR)/lj_cl.h;

$(OBJ_DIR)/lj.o: $(ALL_H) lj.h lj.cpp  $(OBJ_DIR)/lj_cl.h $(OBJ_DIR)/lj_cl.h $(OBJ_DIR)/base_atomic.o
	$(OCL) -o $@ -c lj.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_ext.o: $(ALL_H) lj.h lj_ext.cpp base_atomic.h
	$(OCL) -o $@ -c lj_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_cl.h: lj_coul.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_coul lj_coul.cu $(OBJ_DIR)/lj_coul_cl.h;

$(OBJ_DIR)/lj_coul.o: $(ALL_H) lj_coul.h lj_coul.cpp  $(OBJ_DIR)/lj_coul_cl.h $(OBJ_DIR)/lj_coul_cl.h $(OBJ_DIR)/base_charge.o
	$(OCL) -o $@ -c lj_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_ext.o: $(ALL_H) lj_coul.h lj_coul_ext.cpp base_charge.h
	$(OCL) -o $@ -c lj_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_long_cl.h: lj_coul_long.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_coul_long lj_coul_long.cu $(OBJ_DIR)/lj_coul_long_cl.h;

$(OBJ_DIR)/lj_coul_long.o: $(ALL_H) lj_coul_long.h lj_coul_long.cpp  $(OBJ_DIR)/lj_coul_long_cl.h $(OBJ_DIR)/base_charge.o
	$(OCL) -o $@ -c lj_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_long_ext.o: $(ALL_H) lj_coul_long.h lj_coul_long_ext.cpp base_charge.h
	$(OCL) -o $@ -c lj_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long_cl.h: lj_class2_long.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_class2_long lj_class2_long.cu $(OBJ_DIR)/lj_class2_long_cl.h;

$(OBJ_DIR)/lj_class2_long.o: $(ALL_H) lj_class2_long.h lj_class2_long.cpp  $(OBJ_DIR)/lj_class2_long_cl.h $(OBJ_DIR)/base_charge.o
	$(OCL) -o $@ -c lj_class2_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long_ext.o: $(ALL_H) lj_class2_long.h lj_class2_long_ext.cpp base_charge.h
	$(OCL) -o $@ -c lj_class2_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_long_cl.h: coul_long.cu
	$(BSH) ./geryon/file_to_cstr.sh coul_long coul_long.cu $(OBJ_DIR)/coul_long_cl.h;

$(OBJ_DIR)/coul_long.o: $(ALL_H) coul_long.h coul_long.cpp  $(OBJ_DIR)/coul_long_cl.h $(OBJ_DIR)/base_charge.o
	$(OCL) -o $@ -c coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_long_ext.o: $(ALL_H) coul_long.h coul_long_ext.cpp base_charge.h
	$(OCL) -o $@ -c coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse_cl.h: morse.cu
	$(BSH) ./geryon/file_to_cstr.sh morse morse.cu $(OBJ_DIR)/morse_cl.h;

$(OBJ_DIR)/morse.o: $(ALL_H) morse.h morse.cpp  $(OBJ_DIR)/morse_cl.h $(OBJ_DIR)/morse_cl.h $(OBJ_DIR)/base_atomic.o
	$(OCL) -o $@ -c morse.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse_ext.o: $(ALL_H) morse.h morse_ext.cpp base_atomic.h
	$(OCL) -o $@ -c morse_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/charmm_long_cl.h: charmm_long.cu
	$(BSH) ./geryon/file_to_cstr.sh charmm_long charmm_long.cu $(OBJ_DIR)/charmm_long_cl.h;

$(OBJ_DIR)/charmm_long.o: $(ALL_H) charmm_long.h charmm_long.cpp  $(OBJ_DIR)/charmm_long_cl.h $(OBJ_DIR)/charmm_long_cl.h $(OBJ_DIR)/base_charge.o
	$(OCL) -o $@ -c charmm_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/charmm_long_ext.o: $(ALL_H) charmm_long.h charmm_long_ext.cpp base_charge.h
	$(OCL) -o $@ -c charmm_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_cl.h: lj96.cu
	$(BSH) ./geryon/file_to_cstr.sh lj96 lj96.cu $(OBJ_DIR)/lj96_cl.h;

$(OBJ_DIR)/lj96.o: $(ALL_H) lj96.h lj96.cpp  $(OBJ_DIR)/lj96_cl.h $(OBJ_DIR)/lj96_cl.h $(OBJ_DIR)/base_atomic.o
	$(OCL) -o $@ -c lj96.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96_ext.o: $(ALL_H) lj96.h lj96_ext.cpp base_atomic.h
	$(OCL) -o $@ -c lj96_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_cl.h: lj_expand.cu
	$(BSH) ./geryon/file_to_cstr.sh lj_expand lj_expand.cu $(OBJ_DIR)/lj_expand_cl.h;

$(OBJ_DIR)/lj_expand.o: $(ALL_H) lj_expand.h lj_expand.cpp  $(OBJ_DIR)/lj_expand_cl.h $(OBJ_DIR)/lj_expand_cl.h $(OBJ_DIR)/base_atomic.o
	$(OCL) -o $@ -c lj_expand.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_ext.o: $(ALL_H) lj_expand.h lj_expand_ext.cpp base_atomic.h
	$(OCL) -o $@ -c lj_expand_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_cl.h: cg_cmm.cu
	$(BSH) ./geryon/file_to_cstr.sh cg_cmm cg_cmm.cu $(OBJ_DIR)/cg_cmm_cl.h;

$(OBJ_DIR)/cg_cmm.o: $(ALL_H) cg_cmm.h cg_cmm.cpp  $(OBJ_DIR)/cg_cmm_cl.h $(OBJ_DIR)/cg_cmm_cl.h $(OBJ_DIR)/base_atomic.o
	$(OCL) -o $@ -c cg_cmm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_ext.o: $(ALL_H) cg_cmm.h cg_cmm_ext.cpp base_atomic.h
	$(OCL) -o $@ -c cg_cmm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_long_cl.h: cg_cmm_long.cu
	$(BSH) ./geryon/file_to_cstr.sh cg_cmm_long cg_cmm_long.cu $(OBJ_DIR)/cg_cmm_long_cl.h;

$(OBJ_DIR)/cg_cmm_long.o: $(ALL_H) cg_cmm_long.h cg_cmm_long.cpp  $(OBJ_DIR)/cg_cmm_long_cl.h $(OBJ_DIR)/cg_cmm_long_cl.h $(OBJ_DIR)/base_atomic.o
	$(OCL) -o $@ -c cg_cmm_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/cg_cmm_long_ext.o: $(ALL_H) cg_cmm_long.h cg_cmm_long_ext.cpp base_charge.h
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

