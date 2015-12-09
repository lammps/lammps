default:
	cd ./01_hello_world; \
	make -j 4
	cd ./01_hello_world_lambda; \
	make -j 4
	cd ./02_simple_reduce; \
	make -j 4
	cd ./02_simple_reduce_lambda; \
	make -j 4
	cd ./03_simple_view; \
	make -j 4
	cd ./03_simple_view_lambda; \
	make -j 4
	cd ./04_simple_memoryspaces; \
	make -j 4
	cd ./05_simple_atomics; \
	make -j 4
	cd ./Advanced_Views; \
	make -j 4
	cd ./Algorithms; \
	make -j 4
	cd ./Hierarchical_Parallelism; \
	make -j 4

openmp:
	cd ./01_hello_world; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./01_hello_world_lambda; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./02_simple_reduce; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./02_simple_reduce_lambda; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./03_simple_view; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./03_simple_view_lambda; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./04_simple_memoryspaces; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./05_simple_atomics; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./Advanced_Views; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./Algorithms; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./Hierarchical_Parallelism; \
	make -j 4 KOKKOS_DEVICES=OpenMP

pthreads:
	cd ./01_hello_world; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./01_hello_world_lambda; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./02_simple_reduce; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./02_simple_reduce_lambda; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./03_simple_view; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./03_simple_view_lambda; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./04_simple_memoryspaces; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./05_simple_atomics; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./Advanced_Views; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./Algorithms; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./Hierarchical_Parallelism; \
	make -j 4 KOKKOS_DEVICES=Pthreads

serial:
	cd ./01_hello_world; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./01_hello_world_lambda; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./02_simple_reduce; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./02_simple_reduce_lambda; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./03_simple_view; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./03_simple_view_lambda; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./04_simple_memoryspaces; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./05_simple_atomics; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./Advanced_Views; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./Algorithms; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./Hierarchical_Parallelism; \
	make -j 4 KOKKOS_DEVICES=Serial

cuda:
	cd ./01_hello_world; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./01_hello_world_lambda; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./02_simple_reduce; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./02_simple_reduce_lambda; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./03_simple_view; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./03_simple_view_lambda; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./04_simple_memoryspaces; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./05_simple_atomics; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./Advanced_Views; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./Algorithms; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./Hierarchical_Parallelism; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial

clean:
	cd ./01_hello_world; \
	make clean
	cd ./01_hello_world_lambda; \
	make clean
	cd ./02_simple_reduce; \
	make clean
	cd ./02_simple_reduce_lambda; \
	make clean
	cd ./03_simple_view; \
	make clean
	cd ./03_simple_view_lambda; \
	make clean
	cd ./04_simple_memoryspaces; \
	make clean
	cd ./05_simple_atomics; \
	make clean
	cd ./Advanced_Views; \
	make clean
	cd ./Algorithms; \
	make clean
	cd ./Hierarchical_Parallelism; \
	make clean

