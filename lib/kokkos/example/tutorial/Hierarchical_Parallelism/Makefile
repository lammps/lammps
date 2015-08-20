default:
	cd ./01_thread_teams; \
	make -j 4
	cd ./01_thread_teams_lambda; \
	make -j 4
	cd ./02_nested_parallel_for; \
	make -j 4
	cd ./03_vectorization; \
	make -j 4
	cd ./04_team_scan; \
	make -j 4

openmp:
	cd ./01_thread_teams; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./01_thread_teams_lambda; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./02_nested_parallel_for; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./03_vectorization; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./04_team_scan; \
	make -j 4 KOKKOS_DEVICES=OpenMP

pthreads:
	cd ./01_thread_teams; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./01_thread_teams_lambda; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./02_nested_parallel_for; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./03_vectorization; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./04_team_scan; \
	make -j 4 KOKKOS_DEVICES=Pthreads

serial:
	cd ./01_thread_teams; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./01_thread_teams_lambda; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./02_nested_parallel_for; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./03_vectorization; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./04_team_scan; \
	make -j 4 KOKKOS_DEVICES=Serial

cuda:
	cd ./01_thread_teams; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./01_thread_teams_lambda; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./02_nested_parallel_for; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./03_vectorization; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./04_team_scan; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial

clean:
	cd ./01_thread_teams; \
	make clean
	cd ./01_thread_teams_lambda; \
	make clean
	cd ./02_nested_parallel_for; \
	make clean
	cd ./03_vectorization; \
	make clean
	cd ./04_team_scan; \
	make clean

