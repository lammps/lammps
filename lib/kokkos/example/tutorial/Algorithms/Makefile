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
	mkdir -p 01_random_numbers
	cd ./01_random_numbers; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Algorithms/01_random_numbers/Makefile ${KOKKOS_SETTINGS}

build-insource:
	cd ./01_random_numbers; \
	$(MAKE) build ${KOKKOS_SETTINGS}
test:
	cd ./01_random_numbers; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Algorithms/01_random_numbers/Makefile ${KOKKOS_SETTINGS}

test-insource:
	cd ./01_random_numbers; \
	$(MAKE) test ${KOKKOS_SETTINGS}
clean:
	cd ./01_random_numbers; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Algorithms/01_random_numbers/Makefile ${KOKKOS_SETTINGS}

clean-insource:
	cd ./01_random_numbers; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
