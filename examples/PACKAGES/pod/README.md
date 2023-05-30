### Compile LAMMPS/POD 

Go to `lammps` directory and build with the POD package:

    cd path/to/lammps
    mkdir build-pod
    cd build-pod
    cmake -C ../cmake/presets/basic.cmake -D PKG_ML-POD=on ../cmake
    cmake --build .

### Fit a POD potential for tantalum

Go to `lammps/examples/PACKAGES/pod/Ta` directory and run 

    lmp -in in.podfit

See the README in `lammps/examples/PACKAGES/pod/Ta` for instructions on how to run MD with the potential.

### Examples for other materials

See [https://github.com/cesmix-mit/pod-examples](https://github.com/cesmix-mit/pod-examples)
