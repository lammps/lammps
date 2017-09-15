
ifndef KOKKOS_PATH
  MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
  KOKKOS_PATH = $(subst Makefile,,$(MAKEFILE_PATH))../..
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
	mkdir -p 01_hello_world
	cd ./01_hello_world; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/01_hello_world/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 01_hello_world_lambda
	cd ./01_hello_world_lambda; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/01_hello_world_lambda/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 02_simple_reduce
	cd ./02_simple_reduce; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 02_simple_reduce_lambda
	cd ./02_simple_reduce_lambda; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce_lambda/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 03_simple_view
	cd ./03_simple_view; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/03_simple_view/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 03_simple_view_lambda
	cd ./03_simple_view_lambda; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/03_simple_view_lambda/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 04_simple_memoryspaces
	cd ./04_simple_memoryspaces; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/04_simple_memoryspaces/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 05_simple_atomics
	cd ./05_simple_atomics; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/05_simple_atomics/Makefile ${KOKKOS_SETTINGS}
	mkdir -p Advanced_Views
	cd ./Advanced_Views; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	mkdir -p Algorithms
	cd ./Algorithms; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Algorithms/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	mkdir -p Hierarchical_Parallelism
	cd ./Hierarchical_Parallelism; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'

build-insource:
	cd ./01_hello_world; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	$(MAKE) build KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	$(MAKE) build KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	$(MAKE) build KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
test:
	cd ./01_hello_world; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/01_hello_world/Makefile ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/01_hello_world_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce/Makefile ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/03_simple_view/Makefile ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/03_simple_view_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/04_simple_memoryspaces/Makefile ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/05_simple_atomics/Makefile ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Algorithms/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'

test-insource:
	cd ./01_hello_world; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	$(MAKE) test KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	$(MAKE) test KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	$(MAKE) test KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
clean:
	cd ./01_hello_world; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/01_hello_world/Makefile ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/01_hello_world_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce/Makefile ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/03_simple_view/Makefile ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/03_simple_view_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/04_simple_memoryspaces/Makefile ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/05_simple_atomics/Makefile ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Algorithms/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'

clean-insource:
	cd ./01_hello_world; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	$(MAKE) clean KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	$(MAKE) clean KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	$(MAKE) clean KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
