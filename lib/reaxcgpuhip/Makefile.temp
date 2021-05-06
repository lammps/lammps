#  Makefile template for Static library. 
# 1. Compile every *.cpp in the folder 
# 2. All obj files under obj folder
# 3. static library .a at lib folder
# 4. run 'make dirmake' before calling 'make'


HIPCC = /opt/rocm/bin/hipcc


HIPCCFLAGS  =  -w -D_SINGLE_DOUBLE -DHAVE_CONFIG_H -DUSE_HIP_DEVICE_SORT  -O3  -fpermissive -D__HIP_PLATFORM_HCC__  -pthread

HIPCCINCS = -I/opt/rocm/include -I/opt/rocm/hcc/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include

HIPCCDEFS = $(HIPCCFLAGS) $(HIPCCINCS)


OUT_FILE_NAME = libreaxcgpu.a

OBJ_DIR=./obj

OUT_DIR=./

# Enumerating of every *.cpp as *.o and using that as dependency.	
# filter list of .c files in a directory.
# FILES =dump_l.c \
#	kter.c \
#
# $(OUT_FILE_NAME): $(patsubst %.c,$(OBJ_DIR)/%.o,$(wildcard $(FILES))) 


# Enumerating of every *.cpp as *.o and using that as dependency
$(OUT_FILE_NAME): $(patsubst %.cu,$(OBJ_DIR)/%.o,$(wildcard *.cu))
	ar -rcs -o $(OUT_DIR)/$@ $^



#Compiling every *.cpp to *.o
$(OBJ_DIR)/%.o: %.cu dirmake
	$(HIPCC) -c $(HIPCCDEFS) -o $@  $<
	
dirmake:
	@mkdir -p $(OUT_DIR)
	@mkdir -p $(OBJ_DIR)
	
clean:
	rm -f $(OBJ_DIR)/*.o $(OUT_DIR)/$(OUT_FILE_NAME) Makefile.bak

rebuild: clean build


 
