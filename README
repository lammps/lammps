# Compile LAMMPS/POD 

  1. go to lammps directory
  2. mkdir build
  3. cd build
  4. cmake -C ../cmake/presets/basic.cmake -D BUILD_SHARED_LIBS=on -D LAMMPS_EXCEPTIONS=on -D PKG_PYTHON=on -D PKG_ML-POD=on ../cmake
  5. cmake --build .

# Run an example to fit a POD potential for Tantalum element

  1. Go to lammps/examples/pod/Ta directory
  2.  ../../../build/lmp -in in.podfit -sc tmp

# Run an example to fit a POD potential for InP compound 

  1. Go to lammps/examples/pod/InP directory
  2.  ../../../build/lmp -in in.podfit -sc tmp

# Run an example to fit a POD potential for GaN compound 

  1. Go to lammps/examples/pod/GaN directory
  2.  ../../../build/lmp -in in.podfit -sc tmp
