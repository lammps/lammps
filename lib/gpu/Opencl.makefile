# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
OCL_H  = $(wildcard ./geryon/ocl*.h) $(UCL_H) lal_preprocessor.h
PRE1_H = lal_preprocessor.h lal_aux_fun1.h
ALL_H  =  $(OCL_H) $(wildcard ./lal_*.h)

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

$(OBJ_DIR)/%_cl.h: lal_%.cu $(PRE1_H)
	$(BSH) ./geryon/file_to_cstr.sh $* $(PRE1_H) $< $@;

# host code compilation

$(OBJ_DIR)/lal_%.o: lal_%.cpp $(KERS)
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

