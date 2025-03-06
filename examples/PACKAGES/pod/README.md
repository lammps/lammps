### Compile LAMMPS/POD 

Go to `lammps` directory and build with the POD package:

    cd path/to/lammps
    mkdir build
    cd build
    cmake -C ../cmake/presets/basic.cmake -D PKG_ML-POD=on ../cmake
    cmake --build .

### Compile LAMMPS/POD with Kokkos 

    cmake -C ../cmake/presets/basic.cmake -C ../cmake/presets/kokkos-cuda.cmake -D PKG_ML-POD=on ../cmake

### Fit a POD potential for Tantalum

Go to `lammps/examples/PACKAGES/pod/Ta` directory and run 

    lmp -in Ta_fit.pod

This creates `Ta_coefficients.pod` for the linear model, which we can use to run MD with

    lmp -in Ta_mdrun.pod

### Fit a POD potential for Indium Phosphide

Go to `lammps/examples/PACKAGES/pod/InP` directory and run 

    lmp -in InP_fit.pod

This creates `InP_coefficients.pod` for the linear model, which we can use to run MD with

    lmp -in InP_mdrun.pod

### Examples for other materials

See [https://github.com/cesmix-mit/pod-examples](https://github.com/cesmix-mit/pod-examples)

