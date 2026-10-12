# IPython and Jupyter Notebooks

This folder contains examples showcasing the usage of the LAMMPS Python
interface and Jupyter notebooks. To use this you will need LAMMPS compiled as
a shared library and the LAMMPS Python package installed.

An extensive guide on how to achieve this is documented in the [LAMMPS manual](https://docs.lammps.org/Python_install.html). There is also a [LAMMPS Python tutorial](https://docs.lammps.org/Howto_python.html).

The following will show one way of creating a Python virtual environment
which has both LAMMPS and its Python package installed:

1. Clone the LAMMPS source code

   ```shell
   $ git clone -b stable https://github.com/lammps/lammps.git
   $ cd lammps
   ```

2. Create a virtual environment for Python (here inside the future build folder)

   ```shell
   $ python3 -m venv build/myenv
   ```

3. Extend `LD_LIBRARY_PATH` (Unix/Linux) or `DYLD_LIBRARY_PATH` (MacOS)

   On Unix/Linux:
   ```shell
   $ echo 'export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH' >> build/myenv/bin/activate
   ```

   On MacOS:
   ```shell
   echo 'export DYLD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$DYLD_LIBRARY_PATH' >> build/myenv/bin/activate
   ```

4. Activate the virtual environment

   ```shell
   $ source build/myenv/bin/activate
   (myenv)$
   ```

5. Configure LAMMPS compilation (CMake)

   ```shell
   (myenv)$ cmake -S cmake -B build -C cmake/presets/basic.cmake \
                  -D BUILD_SHARED_LIBS=on \
                  -D PKG_PYTHON=on \
                  -D CMAKE_INSTALL_PREFIX=$VIRTUAL_ENV
   ```

6. Compile LAMMPS

   ```shell
   (myenv)$ cmake --build build
   ```

7. Install LAMMPS and Python package into virtual environment

   ```shell
   (myenv)$ cmake --build build --target install-python
   ```

8. Install other Python packages into virtual environment

   ```shell
   (myenv)$ pip install jupyter matplotlib pandas mpi4py
   ```

9. Navigate to ipython examples folder

   ```shell
   (myenv)$ cd python/examples/ipython
   ```

10. Launch Jupyter and work inside browser

    ```shell
    (myenv)$ jupyter notebook
    ```
