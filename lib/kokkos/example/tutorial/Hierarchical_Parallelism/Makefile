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
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 01_thread_teams_lambda
	cd ./01_thread_teams_lambda; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams_lambda/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 02_nested_parallel_for
	cd ./02_nested_parallel_for; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/02_nested_parallel_for/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 03_vectorization
	cd ./03_vectorization; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/03_vectorization/Makefile ${KOKKOS_SETTINGS}
	mkdir -p 04_team_scan
	cd ./04_team_scan; \
	$(MAKE) build -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/04_team_scan/Makefile ${KOKKOS_SETTINGS}

build-insource:
	cd ./01_thread_teams; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	$(MAKE) build ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	$(MAKE) build ${KOKKOS_SETTINGS}
test:
	cd ./01_thread_teams; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams/Makefile ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/02_nested_parallel_for/Makefile ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/03_vectorization/Makefile ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	$(MAKE) test -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/04_team_scan/Makefile ${KOKKOS_SETTINGS}

test-insource:
	cd ./01_thread_teams; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	$(MAKE) test ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	$(MAKE) test ${KOKKOS_SETTINGS}
clean:
	cd ./01_thread_teams; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams/Makefile ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/01_thread_teams_lambda/Makefile ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/02_nested_parallel_for/Makefile ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/03_vectorization/Makefile ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	$(MAKE) clean -f ${KOKKOS_PATH}/example/tutorial/Hierarchical_Parallelism/04_team_scan/Makefile ${KOKKOS_SETTINGS}

clean-insource:
	cd ./01_thread_teams; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./01_thread_teams_lambda; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./02_nested_parallel_for; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./03_vectorization; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
	cd ./04_team_scan; \
	$(MAKE) clean ${KOKKOS_SETTINGS}
