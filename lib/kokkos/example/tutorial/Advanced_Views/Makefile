ifndef KOKKOS_PATH
  MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
  KOKKOS_PATH = $(subst Makefile,,$(MAKEFILE_PATH))../../..
endif

ifndef KOKKOS_SETTINGS
  KOKKOS_SETTINGS = "KOKKOS_PATH=${KOKKOS_PATH}"
  ifdef KOKKOS_ARCH
    KOKKOS_SETTINGS += "KOKKOS_ARCH=${KOKKOS_ARCH}"
  endif
  ifdef KOKKOS_DEVICES
    KOKKOS_SETTINGS += "KOKKOS_DEVICES=${KOKKOS_DEVICES}"
  endif
  ifdef KOKKOS_OPTIONS
    KOKKOS_SETTINGS += "KOKKOS_OPTIONS=${KOKKOS_OPTIONS}"
  endif
  ifdef KOKKOS_CUDA_OPTIONS
    KOKKOS_SETTINGS += "KOKKOS_CUDA_OPTIONS=${KOKKOS_CUDA_OPTIONS}"
  endif
endif

build:
	mkdir -p 01_data_layouts
	cd ./01_data_layouts; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/01_data_layouts/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 02_memory_traits
	cd ./02_memory_traits; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/02_memory_traits/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 03_subviews
	cd ./03_subviews; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/03_subviews/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 04_dualviews
	cd ./04_dualviews; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/04_dualviews/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 05_NVIDIA_UVM
	cd ./05_NVIDIA_UVM; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/05_NVIDIA_UVM/Makefile ${KOKKOS_SETTINGS}
	#mkdir -p 06_AtomicViews
	#cd ./06_AtomicViews; \
	#make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/06_AtomicViews/Makefile ${KOKKOS_SETTINGS}
	#mkdir -p 07_Overlapping_DeepCopy
	#cd ./07_Overlapping_DeepCopy; \
	#make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/07_Overlapping_DeepCopy/Makefile ${KOKKOS_SETTINGS}

build-insource:
	cd ./01_data_layouts; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	make build -j 4 ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#make build -j 4 ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#make build -j 4 ${KOKKOS_SETTINGS}
test:
	cd ./01_data_layouts; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/01_data_layouts/Makefile ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/02_memory_traits/Makefile ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/03_subviews/Makefile ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/04_dualviews/Makefile ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/05_NVIDIA_UVM/Makefile ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/06_AtomicViews/Makefile ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/07_Overlapping_DeepCopy/Makefile ${KOKKOS_SETTINGS}

test-insource:
	cd ./01_data_layouts; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	make test -j 4 ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#make test -j 4 ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#make test -j 4 ${KOKKOS_SETTINGS}
clean:
	cd ./01_data_layouts; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/01_data_layouts/Makefile ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/02_memory_traits/Makefile ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/03_subviews/Makefile ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/04_dualviews/Makefile ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/05_NVIDIA_UVM/Makefile ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#make clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/06_AtomicViews/Makefile ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#make clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/07_Overlapping_DeepCopy/Makefile ${KOKKOS_SETTINGS}

clean-insource:
	cd ./01_data_layouts; \
	make clean ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	make clean ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	make clean ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	make clean ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	make clean ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#make clean ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#make clean ${KOKKOS_SETTINGS}
