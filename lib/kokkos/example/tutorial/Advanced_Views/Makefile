default:
	cd ./01_data_layouts; \
	make -j 4
	cd ./02_memory_traits; \
	make -j 4
	cd ./03_subviews; \
	make -j 4
	cd ./04_dualviews; \
	make -j 4
	cd ./05_NVIDIA_UVM; \
	make -j 4
	cd ./06_AtomicViews; \
	make -j 4

openmp:
	cd ./01_data_layouts; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./02_memory_traits; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./03_subviews; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./04_dualviews; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./05_NVIDIA_UVM; \
	make -j 4 KOKKOS_DEVICES=OpenMP
	cd ./06_AtomicViews; \
	make -j 4 KOKKOS_DEVICES=OpenMP

pthreads:
	cd ./01_data_layouts; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./02_memory_traits; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./03_subviews; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./04_dualviews; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./05_NVIDIA_UVM; \
	make -j 4 KOKKOS_DEVICES=Pthreads
	cd ./06_AtomicViews; \
	make -j 4 KOKKOS_DEVICES=Pthreads

serial:
	cd ./01_data_layouts; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./02_memory_traits; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./03_subviews; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./04_dualviews; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./05_NVIDIA_UVM; \
	make -j 4 KOKKOS_DEVICES=Serial
	cd ./06_AtomicViews; \
	make -j 4 KOKKOS_DEVICES=Serial

cuda:
	cd ./01_data_layouts; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./02_memory_traits; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./03_subviews; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./04_dualviews; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./05_NVIDIA_UVM; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial
	cd ./06_AtomicViews; \
	make -j 4 KOKKOS_DEVICES=Cuda,Serial

clean:
	cd ./01_data_layouts; \
	make clean
	cd ./02_memory_traits; \
	make clean
	cd ./03_subviews; \
	make clean
	cd ./04_dualviews; \
	make clean
	cd ./05_NVIDIA_UVM; \
	make clean
	cd ./06_AtomicViews; \
	make clean

