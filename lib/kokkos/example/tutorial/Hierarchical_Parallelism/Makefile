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
	mkdir -p 01_thread_teams
	cd ./01_thread_teams; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 01_thread_teams_lambda
	cd ./01_thread_teams_lambda; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams_lambda/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 02_nested_parallel_for
	cd ./02_nested_parallel_for; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/02_nested_parallel_for/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 03_vectorization
	cd ./03_vectorization; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/03_vectorization/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 04_team_scan
	cd ./04_team_scan; \
	make build -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/04_team_scan/Makefile ${KOKKOS_SETTINGS}

build-insource:
	cd ./01_thread_teams; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	make build -j 4 ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	make build -j 4 ${KOKKOS_SETTINGS}
test:
	cd ./01_thread_teams; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams/Makefile ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/02_nested_parallel_for/Makefile ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/03_vectorization/Makefile ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	make test -j 4 -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/04_team_scan/Makefile ${KOKKOS_SETTINGS}

test-insource:
	cd ./01_thread_teams; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	make test -j 4 ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	make test -j 4 ${KOKKOS_SETTINGS}
clean:
	cd ./01_thread_teams; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams/Makefile ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/02_nested_parallel_for/Makefile ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/03_vectorization/Makefile ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	make clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/04_team_scan/Makefile ${KOKKOS_SETTINGS}

clean-insource:
	cd ./01_thread_teams; \
	make clean ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	make clean ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	make clean ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	make clean ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	make clean ${KOKKOS_SETTINGS}
