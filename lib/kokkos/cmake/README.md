![Kokkos](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# Developing Kokkos

This document contains a build system overview for developers with information on adding new CMake options that could influence
* Header configuration macros
* Optional features
* Third-partly libraries
* Compiler and linker flags
For build system details for users, refer to the [build instructions](../BUILD.md).

## Build System

Kokkos uses CMake to configure, build, and install.
Rather than being a completely straightforward use of modern CMake,
Kokkos has several extra complications, primarily due to:
* Kokkos must support linking to an installed version or in-tree builds as a subdirectory of a larger project.
* Kokkos must configure a special compiler `nvcc_wrapper` that allows `nvcc` to accept all C++ flags (which `nvcc` currently does not).
* Kokkos must work as a part of TriBITS, a CMake library providing a particular build idiom for Trilinos.
* Kokkos has many pre-existing users. We need to be careful about breaking previous versions or generating meaningful error messags if we do break backwards compatibility.

If you are looking at the build system code wondering why certain decisions were made: we have had to balance many competing requirements and certain technical debt. Everything in the build system was done for a reason, trying to adhere as closely as possible to modern CMake best practices while meeting all pre-existing. customer requirements.

### Modern CMake Philosophy

Modern CMake relies on understanding the principle of *building* and *using* a code project.
What preprocessor, compiler, and linker flags do I need to *build* my project?
What flags does a downstream project that links to me need to *use* my project?
In CMake terms, flags that are only needed for building are `PRIVATE`.
Only Kokkos needs these flags, not a package that depends on Kokkos.
Flags that must be used in a downstream project are `PUBLIC`.
Kokkos must tell other projects to use them.

In Kokkos, almost everything is a public flag since Kokkos is driven by headers and Kokkos is in charge of optimizing your code to achieve performance portability!
Include paths, C++ standard flags, architecture-specific optimizations, or OpenMP and CUDA flags are all examples of flags that Kokkos configures and adds to your project.

Modern CMake now automatically propagates flags through the `target_link_libraries` command.
Suppose you have a library `stencil` that needs to build with Kokkos.
Consider the following CMake code:

````
find_package(Kokkos)
add_library(stencil stencil.cpp)
target_link_libraries(stencil Kokkos::kokkos)
````

This locates the Kokkos package, adds your library, and tells CMake to link Kokkos to your library.
All public build flags get added automatically through the `target_link_libraries` command.
There is nothing to do. You can be happily oblivious to how Kokkos was configured.
Everything should just work.

As a Kokkos developer who wants to add new public compiler flags, how do you ensure that CMake does this properly? Modern CMake works through targets and properties.
Each target has a set of standard properties:
* `INTERFACE_COMPILE_OPTIONS` contains all the compiler options that Kokkos should add to downstream projects
* `INTERFACE_INCLUDE_DIRECTORIES` contains all the directories downstream projects must include from Kokkos
* `INTERFACE_COMPILE_DEFINITIONS` contains the list of preprocessor `-D` flags
* `INTERFACE_LINK_LIBRARIES` contains all the libraries downstream projects need to link
* `INTERFACE_COMPILE_FEATURES` essentially adds compiler flags, but with extra complications. Features names are specific to CMake. More later.

CMake makes it easy to append to these properties using:
* `target_compile_options(kokkos PUBLIC -fmyflag)`
* `target_include_directories(kokkos PUBLIC mySpecialFolder)`
* `target_compile_definitions(kokkos PUBLIC -DmySpecialFlag=0)`
* `target_link_libraries(kokkos PUBLIC mySpecialLibrary)`
* `target_compile_features(kokkos PUBLIC mySpecialFeature)`
Note that all of these use `PUBLIC`! Almost every Kokkos flag is not private to Kokkos, but must also be used by downstream projects.


### Compiler Features and Compiler Options
Compiler options are flags like `-fopenmp` that do not need to be "resolved." 
The flag is either on or off.
Compiler features are more fine-grained and require conflicting requests to be resolved.
Suppose I have
````
add_library(A a.cpp)
target_compile_features(A PUBLIC cxx_std_11)
````
then another target
````
add_library(B b.cpp)
target_compile_features(B PUBLIC cxx_std_14)
target_link_libraries(A B)
````
I have requested two diferent features.
CMake understands the requests and knows that `cxx_std_11` is a subset of `cxx_std_14`.
CMake then picks C++14 for library `B`.
CMake would not have been able to do feature resolution if we had directly done:
````
target_compile_options(A PUBLIC -std=c++11)
````

### Adding Kokkos Options
After configuring for the first time,
CMake creates a cache of configure variables in `CMakeCache.txt`.
Reconfiguring in the folder "restarts" from those variables.
All flags passed as `-DKokkos_SOME_OPTION=X` to `cmake` become variables in the cache.
All Kokkos options begin with camel case `Kokkos_` followed by an upper case option name.

CMake best practice is to avoid cache variables, if possible.
In essence, you want the minimal amount of state cached between configurations.
And never, ever have behavior influenced by multiple cache variables.
If you want to change the Kokkos configuration, have a single unique variable that needs to be changed.
Never require two cache variables to be changed.

Kokkos provides a function `KOKKOS_OPTION` for defining valid cache-level variables,
proofreading them, and defining local project variables.
The most common variables are called `Kokkos_ENABLE_X`,
for which a helper function `KOKKOS_ENABLE_OPTION` is provided, e.g.
````
KOKKOS_ENABLE_OPTION(TESTS OFF  "Whether to build tests")
````
The function checks if `-DKokkos_ENABLE_TESTS` was given,
whether it was given with the wrong case, e.g. `-DKokkos_Enable_Tests`,
and then defines a regular (non-cache) variable `KOKKOS_ENABLE_TESTS` to `ON` or `OFF`
depending on the given default and whether the option was specified.

### Defining Kokkos Config Macros

Sometimes you may want to add `#define Kokkos_X` macros to the config header.
This is straightforward with CMake.
Suppose you want to define an optional macro `KOKKOS_SUPER_SCIENCE`.
Simply go into `KokkosCore_config.h.in` and add
````
#cmakedefine KOKKOS_SUPER_SCIENCE
````
I can either add
````
KOKKOS_OPTION(SUPER_SCIENCE ON "Whether to do some super science")
````
to directly set the variable as a command-line `-D` option.
Alternatively, based on other logic, I could add to a `CMakeLists.txt`
````
SET(KOKKOS_SUPER_SCIENCE ON)
````
If not set as a command-line option (cache variable), you must make sure the variable is visible in the top-level scope.
If set in a function, you would need:
````
SET(KOKKOS_SUPER_SCIENCE ON PARENT_SCOPE)
````

### Third-Party Libraries
In much the same way that compiler flags transitively propagate to dependent projects,
modern CMake allows us to propagate dependent libraries.
If Kokkos depends on, e.g. `hwloc` the downstream project will also need to link `hwloc`.
There are three stages in adding a new third-party library (TPL):
* Finding: find the desired library on the system and verify the installation is correct
* Importing: create a CMake target, if necessary, that is compatible with `target_link_libraries`. This is mostly relevant for TPLs not installed with CMake.
* Exporting: make the desired library visible to downstream projects 

TPLs are somewhat complicated by whether the library was installed with CMake or some other build system.
If CMake, our lives are greatly simplified. We simply use `find_package` to locate the installed CMake project then call `target_link_libraries(kokkoscore PUBLIC/PRIVATE TPL)`. For libaries not installed with CMake, the process is a bit more complex.
It is up to the Kokkos developers to "convert" the library into a CMake target as if it had been installed as a valid modern CMake target with properties.  
There are helper functions for simplifying the process of importing TPLs in Kokkos, but we walk through the process in detail to clearly illustrate the steps involved.

#### TPL Search Order

There are several options for where CMake could try to find a TPL.
If there are multiple installations of the same TPL on the system,
the search order is critical for making sure the correct TPL is found.
There are 3 possibilities that could be used:

1. Default system paths like /usr
1. User-provided paths through options `<NAME>_ROOT` and `Kokkos_<NAME>_DIR`
1. Additional paths not in the CMake default list or provided by the user that Kokkos decides to add. For example, Kokkos may query `nvcc` or `LD_LIBRARY_PATH` for where to find CUDA libraries.

The following is the search order that Kokkos follows. Note: This differs from the default search order used by CMake `find_library` and `find_header`. CMake prefers default system paths over user-provided paths.
For Kokkos (and package managers in general), it is better to prefer user-provided paths since this usually indicates a specific version we want.

1. `<NAME>_ROOT`
1. `Kokkos_<NAME>_DIR`
1.  Paths added by Kokkos CMake logic
1.  Default system paths (if allowed)

Default system paths are allowed in two cases. First, none of the other options are given so the only place to look is system paths. Second, if explicitly given permission, configure will look in system paths.
The rationale for this logic is that if you specify a custom location, you usually *only* want to look in that location.
If you do not find the TPL where you expect it, you should error out rather than grab another random match.


#### Finding TPLs

If finding a TPL that is not a modern CMake project, refer to the `FindHWLOC.cmake` file in `cmake/Modules` for an example.
You will ususally need to verify expected headers with `find_path`
````
find_path(TPL_INCLUDE_DIR mytpl.h PATHS "${KOKKOS_MYTPL_DIR}/include")
````
This insures that the library header is in the expected include directory and defines the variable `TPL_INCLUDE_DIR` with a valid path if successful.
Similarly, you can verify a library
````
find_library(TPL_LIBRARY mytpl PATHS "${KOKKOS_MYTPL_DIR/lib")
````
that then defines the variable `TPL_LIBRARY` with a valid path if successful.
CMake provides a utility for checking if the `find_path` and `find_library` calls were successful that emulates the behavior of `find_package` for a CMake target.
````
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MYTPL DEFAULT_MSG
                                  MYTPL_INCLUDE_DIR MYTPL_LIBRARY)
````
If the find failed, CMake will print standard error messages explaining the failure.

#### Importing TPLs

The installed TPL must be adapted into a CMake target.
CMake allows libraries to be added that are built externally as follows:
````
add_library(Kokkos::mytpl UNKNOWN IMPORTED)
````
Importantly, we use a `Kokkos::` namespace to avoid name conflicts and identify this specifically as the version imported by Kokkos.
Because we are importing a non-CMake target, we must populate all the target properties that would have been automatically populated for a CMake target.
````
set_target_properties(Kokkos::mytpl PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${MYTPL_INCLUDE_DIR}"
  IMPORTED_LOCATION "${MYTPL_LIBRARY}"
)
````

#### Exporting TPLs

Kokkos may now depend on the target `Kokkos::mytpl` as a `PUBLIC` library (remember building and using).
This means that downstream projects must also know about `Kokkos::myptl` - so Kokkos must export them.
In the `KokkosConfig.cmake.in` file, we need to add code like the following:
````
set(MYTPL_LIBRARY @MYTPL_LIBRARY@)
set(MYTPL_INCLUDE_DIR @MYTPL_INCLUDE_DIR@)
add_library(Kokkos::mytpl UNKNOWN IMPORTED)
set_target_properties(Kokkos::mytpl PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${MYTPL_INCLUDE_DIR}"
  IMPORTED_LOCATION "${MYTPL_LIBRARY}"
)
````
If this looks familiar, that's because it is exactly the same code as above for importing the TPL.
Exporting a TPL really just means importing the TPL when Kokkos is loaded by an external project.
We will describe helper functions that simplify this process.

#### Interface TPLs

If a TPL is just a library and set of headers, we can make a simple `IMPORTED` target.
However, a TPL is actually completely flexible and need not be limited to just headers and libraries.
TPLs can configure compiler flags, linker flags, or multiple different libraries.
For this, we use a special type of CMake target: `INTERFACE` libraries.
These libraries don't build anything.
They simply populate properties that will configure flags for dependent targets.
We consider the example:
````
add_library(PTHREAD INTERFACE)
target_compile_options(PTHREAD PUBLIC -pthread)
````
Kokkos uses the compiler flag `-pthread` to define compiler macros for re-entrant functions rather than treating it simply as a library with header `pthread.h` and library `-lpthread`.
Any property can be configured, e.g.
````
target_link_libraries(MYTPL ...)
````
In contrast to imported TPLs which require direct modification of `KokkosConfig.cmake.in`,
we can use CMake's built-in export functions:
````
INSTALL(
  TARGETS MYTPL
  EXPORT KokkosTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
````
These interface targets will be automatically populated in the config file.

#### Linking the TPL
After finishing the import process, it still remains to link the imported target as needed.
For example,
````
target_link_libraries(kokkoscore PUBLIC Kokkos::HWLOC)
````
The complexity of which includes, options, and libraries the TPL requires
should be encapsulated in the CMake target.

#### TPL Helper Functions
##### KOKKOS_IMPORT_TPL
This function can be invoked as, e.g.
````
KOKKOS_IMPORT_TPL(HWLOC)
````
This function checks if the TPL was enabled by a `-DKokkos_ENABLE_HWLOC=On` flag.
If so, it calls `find_package(TPLHWLOC)`.
This invokes the file `FindTPLHWLOC.cmake` which should be contained in the `cmake/Modules` folder.
If successful, another function `KOKKOS_EXPORT_CMAKE_TPL` gets invoked.
This automatically adds all the necessary import commands to `KokkosConfig.cmake`.

##### KOKKOS_FIND_IMPORTED
Inside a `FindTPLX.cmake` file, the simplest way to import a library is to call, e.g.
````
KOKKOS_FIND_IMPORTED(HWLOC LIBRARY hwloc HEADER hwloc.h)
````
This finds the location of the library and header and creates an imported target `Kokkos::HWLOC`
that can be linked against.
The library/header find can be guided with `-DHWLOC_ROOT=` or `-DKokkos_HWLOC_DIR=` during CMake configure.
These both specify the install prefix.

##### KOKKOS_LINK_TPL
This function checks if the TPL has been enabled.
If so, it links a given library against the imported (or interface) TPL target.

##### KOKKOS_CREATE_IMPORTED_TPL
This helper function is best understood by reading the actual code.
This function takes arguments specifying the properties and creates the actual TPL target.
The most important thing to understand for this function is whether you call this function with the optional `INTERFACE` keyword.
This tells the project to either create the target as an imported target or interface target, as discussed above.

##### KOKKOS_EXPORT_CMAKE_TPL
Even if the TPL just loads a valid CMake target, we still must "export" it into the config file.
When Kokkos is loaded by a downstream project, this TPL must be loaded.
Calling this function simply appends text recording the location where the TPL was found
and adding a `find_dependency(...)` call that will reload the CMake target.

### The Great TriBITS Compromise

TriBITS was a masterpiece of CMake version 2 before the modern CMake idioms of building and using.
TriBITS greatly limited verbosity of CMake files, handled complicated dependency trees between packages, and handled automatically setting up include and linker paths for dependent libraries.

Kokkos is now used by numerous projects that don't (and won't) depend on TriBITS for their build systems.
Kokkos has to work outside of TriBITS and provide a standard CMake 3+ build system.
At the same time, Kokkos is used by numerous projects that depend on TriBITS and don't (and won't) switch to a standard CMake 3+ build system.

Instead of calling functions `TRIBITS_X(...)`, the CMake calls wrapper functions `KOKKOS_X(...)`.
If TriBITS is available (as in Trilinos), `KOKKOS_X` will just be a thin wrapper around `TRIBITS_X`.
If TriBITS is not available, Kokkos maps `KOKKOS_X` calls to native CMake that complies with CMake 3 idioms.
For the time being, this seems the most sensible way to handle the competing requirements of a standalone modern CMake and TriBITS build system.

##### [LICENSE](https://github.com/kokkos/kokkos/blob/devel/LICENSE)

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.
