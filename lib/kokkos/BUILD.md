![Kokkos](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# Installing and Using Kokkos

## Kokkos Philosophy
Kokkos provides a modern CMake style build system.
As C++ continues to develop for C++20 and beyond, CMake is likely to provide the most robust support
for C++.  Applications heavily leveraging Kokkos are strongly encouraged to use a CMake build system.

You can either use Kokkos as an installed package (encouraged) or use Kokkos in-tree in your project.
Modern CMake is exceedingly simple at a high-level (with the devil in the details).
Once Kokkos is installed In your `CMakeLists.txt` simply use:
````cmake
find_package(Kokkos REQUIRED)
````
Then for every executable or library in your project:
````cmake
target_link_libraries(myTarget Kokkos::kokkos)
````
That's it! There is no checking Kokkos preprocessor, compiler, or linker flags.
Kokkos propagates all the necessary flags to your project.
This means not only is linking to Kokkos easy, but Kokkos itself can actually configure compiler and linker flags for *your*
project.
When configuring your project just set:
````bash
> cmake ${srcdir} \
  -DKokkos_ROOT=${kokkos_install_prefix} \
  -DCMAKE_CXX_COMPILER=${compiler_used_to_build_kokkos}
````
Note: You may need the following if your project requires a minimum CMake version older than 3.12:
````cmake
cmake_policy(SET CMP0074 NEW)
````
If building in-tree, there is no `find_package`. You can use `add_subdirectory(kokkos)` with the Kokkos source and again just link with `target_link_libraries(Kokkos::kokkos)`.
The examples in `examples/cmake_build_installed` and `examples/cmake_build_in_tree` can help get you started.


## Configuring CMake
A very basic installation of Kokkos is done with:
````bash
> cmake ${srcdir} \
 -DCMAKE_CXX_COMPILER=g++ \
 -DCMAKE_INSTALL_PREFIX=${kokkos_install_folder}
````
which builds and installed a default Kokkos when you run `make install`.
There are numerous device backends, options, and architecture-specific optimizations that can be configured, e.g.
````bash
> cmake ${srcdir} \
 -DCMAKE_CXX_COMPILER=g++ \
 -DCMAKE_INSTALL_PREFIX=${kokkos_install_folder} \
 -DKokkos_ENABLE_OPENMP=ON
````
which activates the OpenMP backend. All of the options controlling device backends, options, architectures, and third-party libraries (TPLs) are given below.

Kokkos requires as a minimum C++17, however C++20 and C++23 are supported depending on the compiler.

The latest minimum compiler versions can be found in `cmake/kokkos_compiler_id.cmake`.

## Known Issues<a name="KnownIssues"></a>

### Cray

* The Cray compiler wrappers do static linking by default. This seems to break the Kokkos build. You will likely need to set the environment variable `CRAYPE_LINK_TYPE=dynamic` in order to link correctly. Kokkos warns during configure if this is missing.
* The Cray compiler identifies to CMake as Clang, but it sometimes has its own flags that differ from Clang. We try to include all exceptions, but flag errors may occur in which a Clang-specific flag is passed that the Cray compiler does not recognize.

### Fortran

* In a mixed C++/Fortran code, CMake will use the C++ linker by default. If you override this behavior and use Fortran as the link language, the link may break because Kokkos adds linker flags expecting the linker to be C++. Prior to CMake 3.18, Kokkos has no way of detecting in downstream projects that the linker was changed to Fortran.  From CMake 3.18, Kokkos can use generator expressions to avoid adding flags when the linker is not C++. Note: Kokkos will not add any linker flags in this Fortran case. The user will be entirely on their own to add the appropriate linker flags.

## Spack
An alternative to manually building with the CMake is to use the Spack package manager.
Make sure you have downloaded [Spack](https://github.com/spack/spack).
The easiest way to configure the Spack environment is:
````bash
> source spack/share/spack/setup-env.sh
````
with other scripts available for other shells.
You can display information about how to install packages with:
````bash
> spack info kokkos
A basic installation would be done as:
````bash
> spack install kokkos
````
Spack allows options and and compilers to be tuned in the install command.
````bash
> spack install kokkos@3.0 %gcc@7.3.0 +openmp
````
This example illustrates the three most common parameters to Spack:
* Variants: specified with, e.g. `+openmp`, this activates (or deactivates with, e.g. `~openmp`) certain options.
* Version:  immediately following `kokkos` the `@version` can specify a particular Kokkos to build
* Compiler: a default compiler will be chosen if not specified, but an exact compiler version can be given with the `%`option.

For a complete list of Kokkos options, run:
````bash
> spack info kokkos
````
More details can be found in the [Spack README](Spack.md)

#### Spack Development
Spack currently installs packages to a location determined by a unique hash. This hash name is not really "human readable".
Generally, Spack usage should never really require you to reference the computer-generated unique install folder.
If you must know, you can locate Spack Kokkos installations with:
````bash
> spack find -p kokkos ...
````
where `...` is the unique spec identifying the particular Kokkos configuration and version.

A better way to use Spack for doing Kokkos development is the dev-build feature of Spack.
For dev-build details, consult the kokkos-spack repository [README](https://github.com/kokkos/kokkos-spack/blob/master/README.md).

# Kokkos Keyword Listing

Please refer to our [wiki](https://kokkos.github.io/kokkos-core-wiki/keywords.html#cmake-keywords).
