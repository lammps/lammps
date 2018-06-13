cmake-buildsystem
-----------------

To use the cmake build system instead of the make-driven one, do:
```
cmake /path/to/lammps/source/cmake
```
(please note the cmake directory as the very end)

To enable package, e.g. GPU do
```
cmake /path/to/lammps/source/cmake -DPKG_GPU=ON
```

cmake has many many options, do get an overview use the curses-based cmake interface, ccmake:
```
ccmake /path/to/lammps/source/cmake
```
(Don't forget to press "g" for generate once you are done with configuring)
