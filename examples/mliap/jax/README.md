# Running JAX from LAMMPS

### Getting started

First make a Python environment with dependencies:

    conda create --name jax python=3.10
    conda activate jax
    # Upgrade pip
    python -m pip install --upgrade pip
    # Install JAX:
    python -m pip install --upgrade "jax[cpu]"
    # Install other dependencies:
    python -m pip install numpy scipy torch scikit-learn virtualenv psutil tabulate mpi4py Cython

Install LAMMPS:

    cd /path/to/lammps
    mkdir build-jax; cd build-jax
    cmake ../cmake -DLAMMPS_EXCEPTIONS=yes \
                   -DBUILD_SHARED_LIBS=yes \
                   -DMLIAP_ENABLE_PYTHON=yes \
                   -DPKG_PYTHON=yes \
                   -DPKG_ML-SNAP=yes \
                   -DPKG_ML-IAP=yes \
                   -DPYTHON_EXECUTABLE:FILEPATH=`which python`
    make -j4
    make install-python

### Kokkos install

Use same Python dependencies as above, with some extra changes:

1. Make sure you install cupy properly! E.g. 

        python -m pip install cupy-cuda12x

2. Install JAX for GPU/CUDA:

        python -m pip install --trusted-host storage.googleapis.com --upgrade "jax[cuda12_local]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

3. Install cudNN: https://developer.nvidia.com/cudnn

Install LAMMPS. Take care to change `Kokkos_ARCH_*` flag:

    cmake ../cmake -DLAMMPS_EXCEPTIONS=yes \
                  -DBUILD_SHARED_LIBS=yes \
                  -DPKG_PYTHON=yes \
                  -DPKG_ML-SNAP=yes \
                  -DPKG_ML-IAP=yes \
                  -DMLIAP_ENABLE_PYTHON=yes \
                  -DPKG_KOKKOS=yes \
                  -DKokkos_ARCH_TURING75=yes \
                  -DKokkos_ENABLE_CUDA=yes \
                  -DKokkos_ENABLE_OPENMP=yes \
                  -DCMAKE_CXX_COMPILER=${HOME}/lammps/lib/kokkos/bin/nvcc_wrapper \
                  -DPYTHON_EXECUTABLE:FILEPATH=`which python`
    make -j
    make install-python

Run example:

    mpirun -np 1 lmp -k on g 1 -sf kk -pk kokkos newton on -in in.run

### Deploying JAX models on CPU

Use `deploy_script.py`, which will wrap model with `write_unified_jax`.

    python deploy_script.py

This creates `.pkl` file to be loaded by LAMMPS ML-IAP Unified.

Run LAMMPS with the model:

    mpirun -np P lmp -in in.run

### Deploying JAX models in Kokkos

Use `deploy_script_kokkos.py`, which will wrap model with `write_unified_jax_kokkos`.

    python deploy_script_kokkos.py

This creates `.pkl` file to be loaded by LAMMPS ML-IAP Unified.

Run LAMMPS with the model:

    mpirun -np 1 lmp -k on g 1 -sf kk -pk kokkos newton on -in in.run