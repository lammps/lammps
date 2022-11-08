### Compile LAMMPS/POD 

Go to `lammps` directory and build with the POD package:

    cd path/to/lammps
    mkdir build-pod
    cd build-pod
    cmake -C ../cmake/presets/basic.cmake -D PKG_ML-POD=on ../cmake
    cmake --build .

### Run an example to fit a POD potential for Tantalum element

Go to lammps/examples/pod/Ta directory and run 

    ../../../build-pod/lmp -in in.podfit

Also see the README in the `Ta` directory for instructions on how to run MD with the potential.

### Examples for other materials

See [https://github.com/cesmix-mit/pod-examples](https://github.com/cesmix-mit/pod-examples)
