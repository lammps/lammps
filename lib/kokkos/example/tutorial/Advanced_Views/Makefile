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
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/01_data_layouts/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 02_memory_traits
	cd ./02_memory_traits; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/02_memory_traits/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 03_subviews
	cd ./03_subviews; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/03_subviews/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 04_dualviews
	cd ./04_dualviews; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/04_dualviews/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 05_NVIDIA_UVM
	cd ./05_NVIDIA_UVM; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/05_NVIDIA_UVM/Makefile ${KOKKOS_SETTINGS}
	#mkdir -p 06_AtomicViews
	#cd ./06_AtomicViews; \
	#$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/06_AtomicViews/Makefile ${KOKKOS_SETTINGS}
	#mkdir -p 07_Overlapping_DeepCopy
	#cd ./07_Overlapping_DeepCopy; \
	#$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/07_Overlapping_DeepCopy/Makefile ${KOKKOS_SETTINGS}

build-insource:
	cd ./01_data_layouts; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#$(MAKE) build ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#$(MAKE) build ${KOKKOS_SETTINGS}

test:
	cd ./01_data_layouts; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/01_data_layouts/Makefile ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/02_memory_traits/Makefile ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/03_subviews/Makefile ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/04_dualviews/Makefile ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/05_NVIDIA_UVM/Makefile ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/06_AtomicViews/Makefile ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/07_Overlapping_DeepCopy/Makefile ${KOKKOS_SETTINGS}

test-insource:
	cd ./01_data_layouts; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#$(MAKE) test ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#$(MAKE) test ${KOKKOS_SETTINGS}

clean:
	cd ./01_data_layouts; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/01_data_layouts/Makefile ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/02_memory_traits/Makefile ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/03_subviews/Makefile ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/04_dualviews/Makefile ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/05_NVIDIA_UVM/Makefile ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/06_AtomicViews/Makefile ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/07_Overlapping_DeepCopy/Makefile ${KOKKOS_SETTINGS}

clean-insource:
	cd ./01_data_layouts; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./02_memory_traits; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./03_subviews; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./04_dualviews; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./05_NVIDIA_UVM; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	#cd ./06_AtomicViews; \
	#$(MAKE) clean ${KOKKOS_SETTINGS}
	#cd ./07_Overlapping_DeepCopy; \
	#$(MAKE) clean ${KOKKOS_SETTINGS}
