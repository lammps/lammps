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

Make sure you install cupy properly!

    python -m pip cupy

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

### Wrapping JAX code

Take inspiration from the `FitSNAP` ML-IAP wrapper: https://github.com/rohskopf/FitSNAP/blob/mliap-unified/fitsnap3lib/tools/write_unified.py

First define JAX model in `deploy_script.py`, which will wrap model with `write_unified`.

    python deploy_script.py

This creates `.pkl` file to be loaded by LAMMPS ML-IAP Unified.

Run LAMMPS with the model:

    mpirun -np P lmp -in in.run