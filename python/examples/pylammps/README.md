# PyLammps and Jupyter Notebooks

This folder contains examples showcasing the usage of the PyLammps Python
interface and Jupyter notebooks. To use this you will need LAMMPS compiled as
a shared library and the LAMMPS Python package installed.

An extensive guide on how to achieve this is documented in the [LAMMPS manual](https://docs.lammps.org/Python_install.html). There is also a [PyLammps tutorial](https://docs.lammps.org/Howto_pylammps.html).

The following will show one way of creating a Python virtual environment
which has both LAMMPS and its Python package installed:

1. Clone the LAMMPS source code

   ```shell
   $ git clone -b stable https://github.com/lammps/lammps.git
   $ cd lammps
   ```

2. Create a build folder

   ```shell
   $ mkdir build
   $ cd build
   ```

3. Create a virtual environment for Python

   ```shell
   $ python3 -m venv myenv
   ```

4. Extend `LD_LIBRARY_PATH` (Unix/Linux) or `DYLD_LIBRARY_PATH` (MacOS)

   On Unix/Linux:
   ```shell
   $ echo 'export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH' >> myenv/bin/activate
   ```

   On MacOS:
   ```shell
   echo 'export DYLD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$DYLD_LIBRARY_PATH' >> myenv/bin/activate
   ```

5. Activate the virtual environment

   ```shell
   $ source myenv/bin/activate
   (myenv)$
   ```

6. Configure LAMMPS compilation (CMake)

   ```shell
   (myenv)$ cmake -C ../cmake/presets/basic.cmake \
                  -D BUILD_SHARED_LIBS=on \
                  -D LAMMPS_EXCEPTIONS=on -D PKG_PYTHON=on \
                  -D CMAKE_INSTALL_PREFIX=$VIRTUAL_ENV \
                  ../cmake
   ```

7. Compile LAMMPS

   ```shell
   (myenv)$ cmake --build .
   ```

8. Install LAMMPS and Python package into virtual environment

   ```shell
   (myenv)$ cmake --install .
   ```

9. Install other Python packages into virtual environment

   ```shell
   (myenv)$ pip install jupyter matplotlib mpi4py
   ```

10. Navigate to pylammps examples folder

    ```shell
    (myenv)$ cd ../python/examples/pylammmps
    ```

11. Launch Jupyter and work inside browser

    ```shell
    (myenv)$ jupyter notebook
    ```
