# Compile LAMMPS with MLIAP

Make sure you have a python environment active with JAX installed
pip install Cython numpy cupy

```
cd lammps # root directory

mkdir build_mliap
cd build_mliap
cp ../cmake/presets/kokkos-cuda.cmake ./my_kokkos-cuda.cmake
# Edit my_kokkos-cuda.cmake at least changing the Kokkos_ARCH_ to your own (see https://docs.lammps.org/Build_extras.html#kokkos-package)

cmake -C my_kokkos-cuda.cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=$(pwd) \
  -D PKG_ML-IAP=ON -D MLIAP_ENABLE_PYTHON=ON \
  -D PKG_ML-SNAP=ON -D PKG_PYTHON=ON -D BUILD_SHARED_LIBS=ON \
  ../cmake

# If it succeeds, compile
make -j 8

# To install the package in your env
make install-python
```

# Execute

```
cd python/examples/mliap/jax # in this folder

# install mytest module in edit mode
pip install -e .  

# create a pickle file
python makemodel.py 

../../../../build_mliap/lmp -k on g 1 -sf kk -pk kokkos newton on neigh half -in simple.in
```