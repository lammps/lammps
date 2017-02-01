
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
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/01_hello_world/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 01_hello_world_lambda
	cd ./01_hello_world_lambda; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/01_hello_world_lambda/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 02_simple_reduce
	cd ./02_simple_reduce; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 02_simple_reduce_lambda
	cd ./02_simple_reduce_lambda; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce_lambda/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 03_simple_view
	cd ./03_simple_view; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/03_simple_view/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 03_simple_view_lambda
	cd ./03_simple_view_lambda; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/03_simple_view_lambda/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 04_simple_memoryspaces
	cd ./04_simple_memoryspaces; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/04_simple_memoryspaces/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 05_simple_atomics
	cd ./05_simple_atomics; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/05_simple_atomics/Makefile ${KOKKOS_SETTINGS}
	mkdir -p Advanced_Views
	cd ./Advanced_Views; \
	make build -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	mkdir -p Algorithms
	cd ./Algorithms; \
	make build -f ${KOKKOS_PATH}/example/tutorial/Algorithms/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	mkdir -p Hierarchical_Parallelism
	cd ./Hierarchical_Parallelism; \
	make build -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'

build-insource:
	cd ./01_hello_world; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	make build KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	make build KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	make build KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
test:
	cd ./01_hello_world; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/01_hello_world/Makefile ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/01_hello_world_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce/Makefile ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/03_simple_view/Makefile ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/03_simple_view_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/04_simple_memoryspaces/Makefile ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/05_simple_atomics/Makefile ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	make test -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	make test -f ${KOKKOS_PATH}/example/tutorial/Algorithms/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	make test -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'

test-insource:
	cd ./01_hello_world; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	make test KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	make test KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	make test KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
clean:
	cd ./01_hello_world; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/01_hello_world/Makefile ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/01_hello_world_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce/Makefile ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/02_simple_reduce_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/03_simple_view/Makefile ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/03_simple_view_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/04_simple_memoryspaces/Makefile ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/05_simple_atomics/Makefile ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Advanced_Views/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Algorithms/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'

clean-insource:
	cd ./01_hello_world; \
	make clean ${KOKKOS_SETTINGS}
	cd ./01_hello_world_lambda; \
	make clean ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce; \
	make clean ${KOKKOS_SETTINGS}
	cd ./02_simple_reduce_lambda; \
	make clean ${KOKKOS_SETTINGS}
	cd ./03_simple_view; \
	make clean ${KOKKOS_SETTINGS}
	cd ./03_simple_view_lambda; \
	make clean ${KOKKOS_SETTINGS}
	cd ./04_simple_memoryspaces; \
	make clean ${KOKKOS_SETTINGS}
	cd ./05_simple_atomics; \
	make clean ${KOKKOS_SETTINGS}
	cd ./Advanced_Views; \
	make clean KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Algorithms; \
	make clean KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
	cd ./Hierarchical_Parallelism; \
	make clean KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'
