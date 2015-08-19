default:
	cd ./01_random_numbers; \
	make -j 4

openmp:
	cd ./01_random_numbers; \
	make -j 4 KOKKOS_DEVICES=OpenMP

pthreads:
	cd ./01_random_numbers; \
	make -j 4 KOKKOS_DEVICES=Pthreads

serial:
	cd ./01_random_numbers; \
	make -j 4 KOKKOS_DEVICES=Serial

cuda:
	cd ./01_random_numbers; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial

clean:
	cd ./01_random_numbers; \
	make clean

