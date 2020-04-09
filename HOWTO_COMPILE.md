# Compilation with `cmake`

`mkdir build; cd build`

`cmake ../cmake`

`make -j`

NOTE: If there are compilation error messages regarding C++11 standards, run 
`cmake ../cmake -D CMAKE_CXX_FLAGS="-std=c++11"`

# Compilation with `make`

`cd src`

`make yes-user-pace`

`make serial -j`

NOTE: If there are compilation error messages regarding C++11 standards, 
extend the following line in `src/MAKE/Makefile.serial` or `src/MAKE/Makefile.mpi`:

`CCFLAGS =	-g -O3 -std=c++11`

# Updating the pace implementation

`cp <ACE-EVALUATOR-REPO-ROOT>/src/*  <LAMMPS-REPO-ROOT>/src/USER-PACE`