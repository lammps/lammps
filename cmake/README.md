# Building LAMMPS using CMake

LAMMPS recently acquired support for building with CMake thanks to the efforts
of Christoph Junghans (LANL) and Richard Berger (Temple U). One of the key
strengths of CMake is that it can generate the necessary build system files of
your own personal preference. It also enables using common development IDEs such
as Eclipse, Visual Studio, QtCreator, Xcode and many more for LAMMPS
development.

CMake can both be used as a command-line (CLI) utility `cmake` or through one of
its GUIs. `ccmake` is a text-based ui to configure and build CMake projects.
`cmake-gui` is a graphical frontend for running CMake operations that works on
Linux desktop environments, Windows and MacOS X.

The following is a tutorial-style introduction in using the CMake system. It
should give you the necessary foundation to understand how to do the most common
tasks, act as a reference and provide examples of typical use cases.

## Table of Contents

   * [Quick Start for the Impatient](#quick-start-for-the-impatient)
   * [Building LAMMPS using cmake](#building-lammps-using-cmake-1)
      * [Prerequisites](#prerequisites)
      * [Build directory vs. Source Directory](#build-directory-vs-source-directory)
   * [Defining and using presets](#defining-and-using-presets)
   * [Reference](#reference)
      * [Common CMake Configuration Options](#common-cmake-configuration-options)
      * [LAMMPS Configuration Options](#lammps-configuration-options)
      * [Parallelization and Accelerator Packages](#parallelization-and-accelerator-packages)
      * [Default Packages](#default-packages)
      * [Other Packages](#other-packages)
      * [User Packages](#user-packages)
      * [Package-Specific Configuration Options](#package-specific-configuration-options)
         * [KSPACE Package](#kspace-package)
         * [MKL](#mkl)
         * [FFTW3](#fftw3)
         * [BLAS](#blas)
         * [LAPACK](#lapack)
         * [PYTHON Package](#python-package)
         * [GPU Package](#gpu-package)
         * [MESSAGE Package](#message-package)
         * [MSCG Package](#mscg-package)
         * [VORONOI Package](#voronoi-package)
         * [USER-LATTE Package](#user-latte-package)
         * [USER-PLUMED Package](#user-plumed-package)
         * [USER-SCAFACOS Package](#user-scafacos-package)
         * [USER-SMD Package](#user-smd-package)
      * [Optional Features](#optional-features)
         * [zlib support](#zlib-support)
         * [JPEG support](#jpeg-support)
         * [PNG support](#png-support)
         * [GZIP support](#gzip-support)
         * [FFMPEG support](#ffmpeg-support)
      * [Compilers](#compilers)
         * [Building with GNU Compilers](#building-with-gnu-compilers)
         * [Building with Intel Compilers](#building-with-intel-compilers)
         * [Building with LLVM/Clang Compilers](#building-with-llvmclang-compilers)

## Quick Start for the Impatient
If you want to skip ahead and just run the compilation using `cmake`, please
find a minimal example below. Together with the options reference below, this
should get you started.

```bash
git clone https://github.com/lammps/lammps.git
mkdir lammps/build
cd lammps/build
cmake [-D OPTION_A=VALUE_A -D OPTION_B=VALUE_B ...] ../cmake
make
```

# Building LAMMPS using `cmake`

## Prerequisites
This tutorial assumes you are running in a command-line environment using a
shell like Bash.

* Linux: any terminal window will work
* MacOS X: launch the Terminal app
* Windows 10: install and run "Bash on Windows" (aka Ubuntu on Windows)

Before we start, please download the latest and greatest version of LAMMPS from
GitHub. You can either download it as a tarball or ZIP file, or via git. When
you start with a fresh lammps directory, the contents should look like this:

```bash
$ ls
bench  doc       lib      potentials  README  tools
cmake  examples  LICENSE  python      src
```

## Build directory vs. Source Directory

By using CMake we separate building LAMMPS into multiple phases:

1. **Configuration**: define which features we want to enable/disable and how it should be compiled
2. **Compilation**: compile each source file and generate binaries and libraries based on the configuration
3. **Installation** (Optional): finally we can install the generated binaries on our system

In the GNU Make based build system of LAMMPS, configuration occurs by running
special make targets like `make yes-MOLECULAR`. These targets modify the
**source directory** (`src/`) directory by copying package files and patching
Makefile. In some cases you are force to manually edit Makefiles to add compiler
options and/or correct include directory and library paths.

These edits and copy operations are no longer necessary when compiling with
CMake. The source directory stays untouched, so you compile LAMMPS in many
different variants using the same source code checkout. It enables true
**out-of-source** builds.

When using Cmake, you can compile in **any** folder outside of the source
directory. Any working directory you choose becomes a so-called **build
directory**. All configuration files and compilation results are stored in this
folder. We recommend calling it something like `build/`.

Let's have a look a quick example, where we get the greatest and latest version
of LAMMPS from GitHub via git:
```bash
git clone https://github.com/lammps/lammps.git
```

We then create a new `build` folder and make it our current working directory:
```
mkdir lammps/build
cd lammps/build
```

To configure LAMMPS we run `cmake` inside of this folder. However it requires at
least one argument. `cmake` needs to read the LAMMPS `CMakeLists.txt` file to
know what to do.  This file is located in the `cmake/` subdirectory of the
lammps checkout. To run `cmake` add the relative or absolute path to the `cmake/`
directory  as first argument.

E.g., if the current working directory is `lammps/build` you can specify the
relative path to `lammps/cmake` as follows:
```
cmake ../cmake
```

You could also specify the absolute path:
```
cmake /path/to/my/lammps/folder/cmake
```

Please note: **This does NOT compile the code!** Running cmake only configures
the next build. It generates the necessary files to compile the code. On
Unix/Linux it defaults to generating Makefiles. You can also choose other output
formats to generate files for Eclipse, Xcode or Visual Studio which are
supported on other platorms.

To compile LAMMPS once the Makefiles are generated, simply type `make` in the
build directory.

```
make
```
# Defining and using presets

The CMake build exposes a lot of different options. In the old build system
some of the package selections were possible by using special make target like
`make yes-std` or `make no-lib`. Achieving a similar result with cmake requires
specifying all options manually. This can quickly become a very long command
line that is hard to handle.  While these could be stored in a simple script
file, there is another way of defining "presets" to compile LAMMPS in a certain
way. Since the cmake build process - contrary to the conventional build system -
includes the compilation of the bundled libraries into the standard build process,
the grouping of those presets is somewhat different.

A preset is a regular CMake script file that can use constructs such as
variables, lists and for-loops to manipulate configuration options and create
an [*initial cache*](https://cmake.org/cmake/help/v3.12/manual/cmake.1.html).
Options must be set with the `CACHE` and `FORCE` flag to ensure they are
considered even during a second cmake run.

Such a file can then be passed to cmake via the `-C` flag. Several examples of
presets can be found in the `cmake/presets` folder.

```bash
# build LAMMPS with all packages enabled which don't use external libraries and enable GPU package
mkdir build
cd build
cmake -C ../cmake/presets/all_on.cmake -C ../cmake/presets/nolib.cmake -D PKG_GPU=on ../cmake
```

# Reference

## Common CMake Configuration Options


<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>CMAKE_INSTALL_PREFIX</code></td>
  <td>Install location where LAMMPS files will be copied to. In the Unix/Linux case with Makefiles this controls what `make install` will do.</td>
  <td>
   Default setting is <code>$HOME/.local</code>.
  </td>
</tr>
<tr>
  <td><code>CMAKE_BUILD_TYPE</code></td>
  <td>Controls if debugging symbols are added to the generated binaries</td>
  <td>
  <dl>
  <dt><code>RelWithDebInfo (default)</code></dt>
  <dt><code>Release</code></dt>
  <dt><code>Debug</code></dt>
  <dt><code>MinSizeRel</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>CMAKE_VERBOSE_MAKEFILE</code></td>
  <td>Enable verbose output from Makefile builds (useful for debugging), the same can be achieved by adding `VERBOSE=1` to the `make` call.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>


## LAMMPS Configuration Options

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>LAMMPS_SIZE_LIMIT</code></td>
  <td>Controls the integer sizes used by LAMMPS internally</td>
  <td>
  <dl>
    <dt><code>LAMMPS_SMALLBIG</code> (default)</dt>
    <dd>32bit , 64bit</dd>
    <dt><code>LAMMPS_SMALLSMALL</code></dt>
    <dd>32bit , 32bit</dd>
    <dt><code>LAMMPS_BIGBIG</code></dt>
    <dd>64bit , 64bit</dd>
  </dl>
  </td>
</tr>
<tr>
  <td><code>LAMMPS_LONGLONG_TO_LONG</code></td>
  <td>Workaround if your system or MPI version does not recognize <code>long long</code> data types</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>LAMMPS_MEMALIGN</code></td>
  <td>controls the alignment of blocks of memory allocated by LAMMPS</td>
  <td>
  <dl>
    <dt><code>64</code> (default)</dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>LAMMPS_EXCEPTIONS</code></td>
  <td>controls whether LAMMPS dies after an error or throws a C++ exception. This is particularly useful when running through the C library interface, since an error
  in LAMMPS then doesn't kill the parent process</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>LAMMPS_MACHINE</code></td>
  <td>allows appending a machine suffix to the generated LAMMPS binary</td>
  <td>
  <dl>
    <dt>*none*  (default)</dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>LAMMPS_LIB_SUFFIX</code></td>
  <td>allows appending a suffix to the generated LAMMPS library</td>
  <td>
  <dl>
    <dt>*none*  (default)</dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>BUILD_LIB</code></td>
  <td>control whether to build LAMMPS as a library</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>BUILD_EXE</code></td>
  <td>control whether to build LAMMPS executable</td>
  <td>
  <dl>
    <dt><code>on</code> (default)</dt>
    <dt><code>off</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>BUILD_SHARED_LIBS</code></td>
  <td>control whether to build LAMMPS as a shared-library</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>BUILD_DOC</code></td>
  <td>control whether to build LAMMPS documentation</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>BUILD_TOOLS</code></td>
  <td>control whether to build LAMMPS tools</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

## Parallelization and Accelerator Packages

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>BUILD_MPI</code></td>
  <td>control whether to build LAMMPS with MPI support. This will look for
  `mpicxx` in your path and use this MPI implementation.</td>
  <td>
  <dl>
    <dt><code>on</code> (default, if found)</dt>
    <dt><code>off</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>BUILD_OMP</code></td>
  <td>control whether to build LAMMPS with OpenMP support.</td>
  <td>
  <dl>
    <dt><code>on</code> (default, if found)</dt>
    <dt><code>off</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_OPT</code></td>
  <td>
  A handful of pair styles which are optimized for improved CPU performance on
  single or multiple cores. These include EAM, LJ, CHARMM, and Morse potentials.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-OMP</code></td>
  <td>
  Hundreds of pair, fix, compute, bond, angle, dihedral, improper, and kspace
  styles which are altered to enable threading on many-core CPUs via OpenMP
  directives.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-INTEL</code></td>
  <td>
  Dozens of pair, fix, bond, angle, dihedral, improper, and kspace styles which
  are optimized for Intel CPUs and KNLs (Knights Landing).
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_GPU</code></td>
  <td>
  Dozens of pair styles and a version of the PPPM long-range Coulombic solver
  optimized for GPUs. All such styles have a “gpu” as a suffix in their style
  name. The GPU code can be compiled with either CUDA or OpenCL, however the
  OpenCL variants are no longer actively maintained and only the CUDA versions
  are regularly tested.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_KOKKOS</code></td>
  <td>Dozens of atom, pair, bond, angle, dihedral, improper, fix, compute styles adapted to compile using the Kokkos library which can convert them to OpenMP or CUDA code so that they run efficiently on multicore CPUs, KNLs, or GPUs.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

## Default Packages

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>PKG_ASPHERE</code></td>
  <td>Computes, time-integration fixes, and pair styles for aspherical particle models including ellipsoids, 2d lines, and 3d triangles.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_BODY</code></td>
  <td>Body-style particles with internal structure. Computes, time-integration fixes, pair styles, as well as the body styles themselves.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_CLASS2</code></td>
  <td>Bond, angle, dihedral, improper, and pair styles for the COMPASS CLASS2 molecular force field.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_COLLOID</code></td>
  <td>Coarse-grained finite-size colloidal particles. Pair styles and fix wall styles for colloidal interactions. Includes the Fast Lubrication Dynamics (FLD) method for hydrodynamic interactions, which is a simplified approximation to Stokesian dynamics.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_COMPRESS</code></td>
  <td>Compressed output of dump files via the zlib compression library, using dump styles with a “gz” in their style name.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_CORESHELL</code></td>
  <td>Compute and pair styles that implement the adiabatic core/shell model for polarizability. The pair styles augment Born, Buckingham, and Lennard-Jones styles with core/shell capabilities. The compute temp/cs command calculates the temperature of a system with core/shell particles.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_DIPOLE</code></td>
  <td>An atom style and several pair styles for point dipole models with short-range or long-range interactions.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_GRANULAR</code></td>
  <td>Pair styles and fixes for finite-size granular particles, which interact with each other and boundaries via frictional and dissipative potentials.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_KSPACE</code></td>
  <td>A variety of long-range Coulombic solvers, as well as pair styles which compute the corresponding short-range pairwise Coulombic interactions. These include Ewald, particle-particle particle-mesh (PPPM), and multilevel summation method (MSM) solvers.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_MANYBODY</code></td>
  <td>
  A variety of manybody and bond-order potentials. These include (AI)REBO, BOP,
  EAM, EIM, Stillinger-Weber, and Tersoff potentials.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_MC</code></td>
  <td>
  Several fixes and a pair style that have Monte Carlo (MC) or MC-like
  attributes. These include fixes for creating, breaking, and swapping bonds,
  for performing atomic swaps, and performing grand-canonical MC (GCMC) in
  conjunction with dynamics.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_MISC</code></td>
  <td>
  A variety of compute, fix, pair, dump styles with specialized capabilities that
  don’t align with other packages.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_MOLECULE</code></td>
  <td>
  A large number of atom, pair, bond, angle, dihedral, improper styles that are
  used to model molecular systems with fixed covalent bonds. The pair styles
  include the Dreiding (hydrogen-bonding) and CHARMM force fields, and a TIP4P
  water model.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_PERI</code></td>
  <td>
  An atom style, several pair styles which implement different Peridynamics
  materials models, and several computes which calculate diagnostics.
  Peridynamics is a a particle-based meshless continuum model.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_QEQ</code></td>
  <td>
  Several fixes for performing charge equilibration (QEq) via different
  algorithms. These can be used with pair styles that perform QEq as part of
  their formulation.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_REPLICA</code></td>
  <td>
  A collection of multi-replica methods which can be used when running multiple
  LAMMPS simulations (replicas). See Section 6.5 for an overview of how to run
  multi-replica simulations in LAMMPS. Methods in the package include nudged
  elastic band (NEB), parallel replica dynamics (PRD), temperature accelerated
  dynamics (TAD), parallel tempering, and a verlet/split algorithm for
  performing long-range Coulombics on one set of processors, and the remainder
  of the force field calcalation on another set.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_RIGID</code></td>
  <td>
  Fixes which enforce rigid constraints on collections of atoms or particles.
  This includes SHAKE and RATTLE, as well as varous rigid-body integrators for a
  few large bodies or many small bodies. Also several computes which calculate
  properties of rigid bodies.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_SHOCK</code></td>
  <td>
  Fixes for running impact simulations where a shock-wave passes through a
  material.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_SPIN</code></td>
  <td>Model atomic magnetic spins classically, coupled to atoms moving in the usual manner via MD. Various pair, fix, and compute styles.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_SNAP</code></td>
  <td>
  A pair style for the spectral neighbor analysis potential (SNAP). SNAP is
  methodology for deriving a highly accurate classical potential fit to a large
  archive of quantum mechanical (DFT) data. Also several computes which analyze
  attributes of the potential.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_SRD</code></td>
  <td>
  A pair of fixes which implement the Stochastic Rotation Dynamics (SRD) method
  for coarse-graining of a solvent, typically around large colloidal particles.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

## Other Packages

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>PKG_KIM</code></td>
  <td>A <code>pair_style kim</code> command which is a wrapper on the Knowledge Base for Interatomic Models (KIM) repository of interatomic potentials, enabling any of them to be used in LAMMPS simulations.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_PYTHON</code></td>
  <td>Enable support for Python scripting inside of LAMMPS.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_MESSAGE</code></td>
  <td>Commands to use LAMMPS as either a client or server and couple it to another application.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_MSCG</code></td>
  <td>
  A fix mscg command which can parameterize a Multi-Scale Coarse-Graining (MSCG)
  model using the open-source MS-CG library.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_MPIIO</code></td>
  <td>
  Support for parallel output/input of dump and restart files via the MPIIO library.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_POEMS</code></td>
  <td>
  A fix that wraps the Parallelizable Open source Efficient Multibody Software
  (POEMS) library, which is able to simulate the dynamics of articulated body
  systems. These are systems with multiple rigid bodies (collections of
  particles) whose motion is coupled by connections at hinge points.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_LATTE</code></td>
  <td>
  A fix command which wraps the LATTE DFTB code, so that molecular dynamics can
  be run with LAMMPS using density-functional tight-binding quantum forces
  calculated by LATTE.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_VORONOI</code></td>
  <td>
  A compute command which calculates the Voronoi tesselation of a collection of atoms by wrapping the Voro++ library. This can be used to calculate the local volume or each atoms or its near neighbors.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

## User Packages

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>PKG_USER-ADIOS</code></td>
  <td>ADIOS is a high-performance I/O library. This package implements the dump “atom/adios” and dump “custom/adios” commands to write data using the ADIOS library.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-ATC</code></td>
  <td>
  ATC stands for atoms-to-continuum. This package implements a fix atc command
  to either couple molecular dynamics with continuum finite element equations or
  perform on-the-fly conversion of atomic information to continuum fields.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-AWPMD</code></td>
  <td>
  AWPMD stands for Antisymmetrized Wave Packet Molecular Dynamics. This package
  implements an atom, pair, and fix style which allows electrons to be treated
  as explicit particles in a classical molecular dynamics model.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-BOCS</code></td>
  <td>
  This package provides fix bocs, a modified version of fix npt which includes the pressure correction to the barostat as outlined in: N. J. H. Dunn and W. G. Noid, “Bottom-up coarse-grained models that accurately describe the structure, pressure, and compressibility of molecular liquids,” J. Chem. Phys. 143, 243148 (2015).
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-CGDNA</code></td>
  <td>
  Several pair styles, a bond style, and integration fixes for coarse-grained
  models of single- and double-stranded DNA based on the oxDNA model of Doye,
  Louis and Ouldridge at the University of Oxford. This includes Langevin-type
  rigid-body integrators with improved stability.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-CGSDK</code></td>
  <td>
  Several pair styles and an angle style which implement the coarse-grained SDK
  model of Shinoda, DeVane, and Klein which enables simulation of ionic liquids,
  electrolytes, lipids and charged amino acids.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-COLVARS</code></td>
  <td>
  COLVARS stands for collective variables, which can be used to implement
  various enhanced sampling methods, including Adaptive Biasing Force,
  Metadynamics, Steered MD, Umbrella Sampling and Restraints. A fix colvars
  command is implemented which wraps a COLVARS library, which implements these
  methods. simulations.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-DIFFRACTION</code></td>
  <td>
  Two computes and a fix for calculating x-ray and electron diffraction
  intensities based on kinematic diffraction theory.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-DPD</code></td>
  <td>
  DPD stands for dissipative particle dynamics. This package implements
  coarse-grained DPD-based models for energetic, reactive molecular crystalline
  materials. It includes many pair styles specific to these systems, including
  for reactive DPD, where each particle has internal state for multiple species
  and a coupled set of chemical reaction ODEs are integrated each timestep.
  Highly accurate time integrators for isothermal, isoenergetic, isobaric and
  isenthalpic conditions are included. These enable long timesteps via the
  Shardlow splitting algorithm.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-DRUDE</code></td>
  <td>
  Fixes, pair styles, and a compute to simulate thermalized Drude oscillators as
  a model of polarization.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-EFF</code></td>
  <td>
  EFF stands for electron force field which allows a classical MD code to model
  electrons as particles of variable radius. This package contains atom, pair,
  fix and compute styles which implement the eFF as described in A.
  Jaramillo-Botero, J. Su, Q. An, and W.A. Goddard III, JCC, 2010. The eFF
  potential was first introduced by Su and Goddard, in 2007.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-FEP</code></td>
  <td>
  FEP stands for free energy perturbation. This package provides methods for
  performing FEP simulations by using a fix adapt/fep command with soft-core
  pair potentials, which have a “soft” in their style name.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-H5MD</code></td>
  <td>
  H5MD stands for HDF5 for MD. HDF5 is a portable, binary, self-describing file
  format, used by many scientific simulations. H5MD is a format for molecular
  simulations, built on top of HDF5. This package implements a dump h5md command
  to output LAMMPS snapshots in this format.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-LB</code></td>
  <td>
  Fixes which implement a background Lattice-Boltzmann (LB) fluid, which can be
  used to model MD particles influenced by hydrodynamic forces.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-MANIFOLD</code></td>
  <td>
  Several fixes and a “manifold” class which enable simulations of particles
  constrained to a manifold (a 2D surface within the 3D simulation box). This is
  done by applying the RATTLE constraint algorithm to formulate single-particle
  constraint functions g(xi,yi,zi) = 0 and their derivative (i.e. the normal of
  the manifold) n = grad(g).
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-MEAMC</code></td>
  <td>
  A pair style for the modified embedded atom (MEAM) potential translated from
  the Fortran version in the MEAM package to plain C++. In contrast to the MEAM
  package, no library needs to be compiled and the pair style can be
  instantiated multiple times.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-MESO</code></td>
  <td>
  Several extensions of the the dissipative particle dynamics (DPD) method.
  Specifically, energy-conserving DPD (eDPD) that can model non-isothermal
  processes, many-body DPD (mDPD) for simulating vapor-liquid coexistence, and
  transport DPD (tDPD) for modeling advection-diffusion-reaction systems. The
  equations of motion of these DPD extensions are integrated through a modified
  velocity-Verlet (MVV) algorithm.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-MGPT</code></td>
  <td>
  A pair style which provides a fast implementation of the quantum-based MGPT
  multi-ion potentials. The MGPT or model GPT method derives from
  first-principles DFT-based generalized pseudopotential theory (GPT) through a
  series of systematic approximations valid for mid-period transition metals
  with nearly half-filled d bands. The MGPT method was originally developed by
  John Moriarty at LLNL. The pair style in this package calculates forces and
  energies using an optimized matrix-MGPT algorithm due to Tomas Oppelstrup at
  LLNL.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-MISC</code></td>
  <td>
  A potpourri of (mostly) unrelated features contributed to LAMMPS by users.
  Each feature is a single fix, compute, pair, bond, angle, dihedral, improper,
  or command style.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-MOFFF</code></td>
  <td>
  Pair, angle and improper styles needed to employ the MOF-FF force field by
  Schmid and coworkers with LAMMPS. MOF-FF is a first principles derived force
  field with the primary aim to simulate MOFs and related porous framework
  materials, using spherical Gaussian charges. It is described in S. Bureekaew
  et al., Phys. Stat. Sol. B 2013, 250, 1128-1141. For the usage of MOF-FF see
  the example in the example directory as well as the MOF+ website.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-MOLFILE</code></td>
  <td>
  A dump molfile command which uses molfile plugins that are bundled with the
  VMD molecular visualization and analysis program, to enable LAMMPS to dump
  snapshots in formats compatible with various molecular simulation tools.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-NETCDF</code></td>
  <td>
  Dump styles for writing NetCDF formatted dump files. NetCDF is a portable,
  binary, self-describing file format developed on top of HDF5. The file
  contents follow the AMBER NetCDF trajectory conventions
  (http://ambermd.org/netcdf/nctraj.xhtml), but include extensions.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-PHONON</code></td>
  <td>
  A fix phonon command that calculates dynamical matrices, which can then be
  used to compute phonon dispersion relations, directly from molecular dynamics
  simulations.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-PLUMED</code></td>
  <td>
  The fix plumed command allows you to use the PLUMED free energy plugin for molecular dynamics to analyze and bias your LAMMPS trajectory on the fly. The PLUMED library is called from within the LAMMPS input script by using the <code>fix plumed</code> command.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-PTM</code></td>
  <td>
  A <code>compute ptm/atom</code> command that calculates local structure characterization using the Polyhedral Template Matching methodology.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-QTB</code></td>
  <td>
  Two fixes which provide a self-consistent quantum treatment of vibrational modes in a classical molecular dynamics simulation. By coupling the MD simulation to a colored thermostat, it introduces zero point energy into the system, altering the energy power spectrum and the heat capacity to account for their quantum nature. This is useful when modeling systems at temperatures lower than their classical limits or when temperatures ramp across the classical limits in a simulation.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-QUIP</code></td>
  <td>
  A pair_style quip command which wraps the QUIP libAtoms library, which
  includes a variety of interatomic potentials, including Gaussian Approximation
  Potential (GAP) models developed by the Cambridge University group.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-QMMM</code></td>
  <td>
  A fix qmmm command which allows LAMMPS to be used in a QM/MM simulation,
  currently only in combination with the Quantum ESPRESSO package.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-REAXC</code></td>
  <td>
  A pair style which implements the ReaxFF potential in C/C++ (in contrast to
  the REAX package and its Fortran library). ReaxFF is universal reactive force
  field. See the src/USER-REAXC/README file for more info on differences between
  the two packages. Also two fixes for monitoring molecules as bonds are created
  and destroyed.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-SCAFACOS</code></td>
  <td>
  A KSpace style which wraps the ScaFaCoS Coulomb solver library to compute long-range Coulombic interactions.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-SDPD</code></td>
  <td>
  A pair style for smoothed dissipative particle dynamics (SDPD), which is an
  extension of smoothed particle hydrodynamics (SPH) to mesoscale where thermal
  fluctuations are important (see the USER-SPH package). Also two fixes for
  moving and rigid body integration of SPH/SDPD particles (particles of
  <code>atom_style meso</code>).</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-SMD</code></td>
  <td>
  An atom style, fixes, computes, and several pair styles which implements
  smoothed Mach dynamics (SMD) for solids, which is a model related to smoothed
  particle hydrodynamics (SPH) for liquids (see the USER-SPH package).
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-SMTBQ</code></td>
  <td>
  A pair style which implements a Second Moment Tight Binding model with QEq
  charge equilibration (SMTBQ) potential for the description of ionocovalent
  bonds in oxides.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-SPH</code></td>
  <td>
  An atom style, fixes, computes, and several pair styles which implements
  smoothed particle hydrodynamics (SPH) for liquids. See the related USER-SMD
  package package for smooth Mach dynamics (SMD) for solids.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-TALLY</code></td>
  <td>
  Several compute styles that can be called when pairwise interactions are
  calculated to tally information (forces, heat flux, energy, stress, etc) about
  individual interactions.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-UEF</code></td>
  <td>
  A fix style for the integration of the equations of motion under extensional
  flow with proper boundary conditions, as well as several supporting compute
  styles and an output option.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-VTK</code></td>
  <td>
  A dump vtk command which outputs snapshot info in the VTK format, enabling
  visualization by Paraview or other visualization packages.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PKG_USER-YAFF</code></td>
  <td>
  Some potentials that are also implemented in the Yet Another Force Field (YAFF) code.
  The expressions and their use are discussed in the following papers:
  <ul>
    <li><a href="https://doi.org/10.1002/jcc.23877" target="_blank">Vanduyfhuys et al., J. Comput. Chem., 36 (13), 1015-1027 (2015)</a></li>
    <li><a href="https://doi.org/10.1002/jcc.25173" target="_blank">Vanduyfhuys et al., J. Comput. Chem., 39 (16), 999-1011 (2018)</a></li>
  </ul>
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

## Package-Specific Configuration Options

### KSPACE Package

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>FFT</code></td>
  <td>
    <p>FFT library for KSPACE package</p>
    <p>If either MKL or FFTW is selected <code>cmake</code> will try to locate
    these libraries automatically. To control which one should be used please see
    the options below for each FFT library.  Otherwise it will default to KISS
    FFT.</p>
  </td>
  <td>
  <dl>
    <dt><code>FFTW3</code></dt>
    <dt><code>FFTW2</code></dt>
    <dt><code>MKL</code></dt>
    <dt><code>KISS</code> (default)</dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>FFT_SINGLE</code></td>
  <td>Use single-precision floating-point in FFT</td>
  <td>
  <dl>
    <dt><code>off</code> (default = double precision)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>FFT_PACK</code></td>
  <td>Optimization for FFT</td>
  <td>
  <dl>
    <dt><code>array (default)</code></dt>
    <dt><code>pointer</code></dt>
    <dt><code>memcpy</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

### FFTW3

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>FFTW3_INCLUDE_DIRS</code></td>
  <td>path to FFTW3 include files</td>
  <td>
  </td>
</tr>
<tr>
  <td><code>FFTW3_LIBRARIES</code></td>
  <td>list of paths to FFTW3 libraries</td>
  <td>
  </td>
</tr>
</tbody>
</table>

### MKL

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>MKL_INCLUDE_DIRS</code></td>
  <td>path to MKL include files</td>
  <td>
  </td>
</tr>
<tr>
  <td><code>MKL_LIBRARIES</code></td>
  <td>list of paths to MKL libraries</td>
  <td>
  </td>
</tr>
</tbody>
</table>

### BLAS

See [FindBLAS documentation](https://cmake.org/cmake/help/latest/module/FindBLAS.html)

### LAPACK

See [FindLAPACK documentation](https://cmake.org/cmake/help/latest/module/FindLAPACK.html)

### PYTHON Package

See [FindPYTHON documentation](https://cmake.org/cmake/help/latest/module/FindPython.html)

### USER-INTEL Package

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>INTEL_ARCH</code></td>
  <td>Target architecture for USER-INTEL package</td>
  <td>
  <dl>
    <dt><code>cpu</code> (default)</dt>
    <dt><code>knl</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>INTEL_LRT_MODE</code></td>
  <td>How to support Long-range thread mode in Verlet integration</td>
  <td>
  <dl>
    <dt><code>threads</code> (default, if pthreads available)</dt>
    <dt><code>none</code> (default, if pthreads not available)</dt>
    <dt><code>c++11</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

### GPU Package
The GPU package builds a support library which can either use OpenCL or CUDA as
target API.

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>GPU_API</code></td>
  <td>API used by GPU package</td>
  <td>
  <dl>
    <dt><code>opencl</code> (default)</dt>
    <dt><code>cuda</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>GPU_PREC</code></td>
  <td>Precision size used by GPU package kernels</td>
  <td>
  <dl>
    <dt><code>mixed</code> (default)</dt>
    <dt><code>single</code></dt>
    <dt><code>double</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>OCL_TUNE</code> (OpenCL only)</td>
  <td>Tuning target for OpenCL driver code</td>
  <td>
  <dl>
    <dt><code>generic</code> (default)</dt>
    <dt><code>intel</code> (Intel CPU)</dt>
    <dt><code>phi</code> (Intel Xeon Phi)</dt>
    <dt><code>fermi</code> (NVIDIA)</dt>
    <dt><code>kepler</code> (NVIDIA)</dt>
    <dt><code>cypress</code> (AMD)</dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>GPU_ARCH</code> (CUDA only)</td>
  <td>CUDA SM architecture targeted by GPU package</td>
  <td>
  <dl>
    <dt><code>sm_20</code> (Fermi)</dt>
    <dt><code>sm_30</code> (Kepler) (default)</dt>
    <dt><code>sm_50</code> (Maxwell)</dt>
    <dt><code>sm_60</code> (Pascal)</dt>
    <dt><code>sm_70</code> (Volta)</dt>
    <dt><code>sm_75</code> (Turing)</dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>CUDPP_OPT</code> (CUDA only)</td>
  <td>Enable CUDA Performance Primitives Optimizations</td>
  <td>
  <dl>
    <dt><code>on</code> (default)</dt>
    <dt><code>off</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>CUDA_MPS_SUPPORT</code> (CUDA only)</td>
  <td>Enable tweaks for running with Nvidia CUDA Multi-process services daemon</td>
  <td>
  <dl>
    <dt><code>on</code></dt>
    <dt><code>off</code> (default)</dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>BIN2C</code> (CUDA only)</td>
  <td>Path to bin2c executable, will automatically pick up the first one in your $PATH.</td>
  <td>(automatic)</td>
</tr>
</tbody>
</table>

### KIM Package

Requires installation of the KIM library with API v2

If `DOWNLOAD_KIM` is set, the KIM library will be downloaded and built inside
the CMake build directory. If the KIM library is already on your system (in a
location CMake cannot find it), set the `PKG_CONFIG_PATH` environment variable
so that `libkim-api` can be found.

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>DOWNLOAD_KIM</code></td>
  <td>Download KIM API v2 and compile it as part of the build.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

### MESSAGE Package

This package can optionally include support for messaging via sockets, using the open-source [ZeroMQ library](http://zeromq.org/), which must be installed on your system.

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>MESSAGE_ZMQ</code></td>
  <td>Build with ZeroMQ support</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ZMQ_LIBRARY</code></td>
  <td>
    ZMQ library file (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
<tr>
  <td><code>ZMG_INCLUDE_DIR</code></td>
  <td>
    Provide include directory of existing ZMQ installation (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
</tbody>
</table>

### MSCG Package

Requires installation of the MSCG library

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>DOWNLOAD_MSCG</code></td>
  <td>Download MSCG and compile it as part of the build</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>MSCG_LIBRARY</code></td>
  <td>
    MSCG library file (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
<tr>
  <td><code>MSCG_INCLUDE_DIR</code></td>
  <td>
    Provide include directory of existing MSCG installation (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
</tbody>
</table>

### VORONOI Package

Requires installation of the Voro++ library

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>DOWNLOAD_VORO</code></td>
  <td>Download Voro++ and compile it as part of the build</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>VORO_LIBRARY</code></td>
  <td>
    Voro++ library file (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
<tr>
  <td><code>VORO_INCLUDE_DIR</code></td>
  <td>
    Provide include directory of existing Voro++ installation (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
</tbody>
</table>

### USER-COLVARS Package

Requires a C++11 compiler to compile with the Lepton library included.

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>COLVARS_LEPTON</code></td>
  <td>Enable the use of the Lepton library inside the Colvars library.
  <td>
  <dl>
    <dt><code>on</code> (default)</dt>
    <dt><code>off</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

### USER-LATTE Package

Requires installation of the LATTE library

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>DOWNLOAD_LATTE</code></td>
  <td>Download LATTE and compile it as part of the build</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>LATTE_LIBRARY</code></td>
  <td>
    LATTE library file (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
</tbody>
</table>

### USER-PLUMED Package

Requires installation of the PLUMED library

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>DOWNLOAD_PLUMED</code></td>
  <td>Download PLUMED and compile it as part of the build</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PLUMED_MODE</code></td>
  <td>
    Determines the linkage mode for the PLUMED library.
  </td>
  <td>
  <dl>
    <dt><code>static</code> (default)</dt>
    <dt><code>shared</code></dt>
    <dt><code>runtime</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>

### USER-LATTE Package

Requires installation of the LATTE library

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>DOWNLOAD_LATTE</code></td>
  <td>Download LATTE and compile it as part of the build</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>LATTE_LIBRARY</code></td>
  <td>
    LATTE library file (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
</tbody>
</table>

### USER-SMD Package

Requires installation of the Eigen3 library

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>DOWNLOAD_EIGEN3</code></td>
  <td>Download Eigen3 and compile it as part of the build</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>EIGEN3_INCLUDE_DIR</code></td>
  <td>
    Provide include directory of existing Eigen3 installation (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
</tbody>
</table>


### USER-SCAFACOS Package

To build with this package, you must download and build the [ScaFaCoS Coulomb solver library](http://www.scafacos.de/)

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>DOWNLOAD_SCAFACOS</code></td>
  <td>Download SCAFACOS and compile it as part of the build</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>SCAFACOS_LIBRARY</code></td>
  <td>
    SCAFACOS library file (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
<tr>
  <td><code>SCAFACOS_INCLUDE_DIR</code></td>
  <td>
    SCAFACOS include directory (only needed if at custom location)
  </td>
  <td>
  </td>
</tr>
</tbody>
</table>

## Optional Features

### zlib support

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>ZLIB_INCLUDE_DIR</code></td>
  <td></td>
  <td>
  </td>
</tr>
<tr>
  <td><code>ZLIB_LIBRARIES</code></td>
  <td></td>
  <td>
  </td>
</tr>
</tbody>
</table>

### JPEG support

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>WITH_JPEG</code></td>
  <td>Enables/Disable JPEG support in LAMMPS</td>
  <td>
  <dl>
    <dt><code>yes</code> (default, if found)</dt>
    <dt><code>no</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>JPEG_INCLUDE_DIR</code></td>
  <td></td>
  <td>
  </td>
</tr>
<tr>
  <td><code>JPEG_LIBRARIES</code></td>
  <td></td>
  <td>
  </td>
</tr>
</tbody>
</table>

### PNG support
(requires zlib support)

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>WITH_PNG</code></td>
  <td>Enables/Disable PNG support in LAMMPS</td>
  <td>
  <dl>
    <dt><code>yes</code> (default, if found)</dt>
    <dt><code>no</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PNG_INCLUDE_DIR</code></td>
  <td></td>
  <td>
  </td>
</tr>
<tr>
  <td><code>PNG_LIBRARIES</code></td>
  <td></td>
  <td>
  </td>
</tr>
</tbody>
</table>

### GZIP support

requires `gzip` to be in your `PATH`

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>WITH_GZIP</code></td>
  <td>Enables/Disable GZIP support in LAMMPS</td>
  <td>
  <dl>
    <dt><code>yes</code> (default, if found)</dt>
    <dt><code>no</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>GZIP_EXECUTABLE</code></td>
  <td>Path to gzip executable, will automatically pick up the first one in your $PATH.</td>
  <td>(automatic)</td>
</tr>
</tbody>
</table>

### FFMPEG support

requires `ffmpeg` to be in your `PATH`

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>WITH_FFMPEG</code></td>
  <td>Enables/Disable FFMPEG support in LAMMPS</td>
  <td>
  <dl>
    <dt><code>yes</code> (default, if found)</dt>
    <dt><code>no</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>FFMPEG_EXECUTABLE</code></td>
  <td>Path to ffmpeg executable, will automatically pick up the first one in your $PATH.</td>
  <td>(automatic)</td>
</tr>
</tbody>
</table>


## Compilers

By default, `cmake` will use your environment C/C++/Fortran compilers for a
build. It uses the `CC`, `CXX` and `FC` environment variables to detect which
compilers should be used. However, these values will be cached after the first
run of `cmake`. Subsequent runs of `cmake` will ignore changes in these
environment variables. To ensure the correct values are used you avoid the
cache by setting the `CMAKE_C_COMPILER`, `CMAKE_CXX_COMPILER`,
`CMAKE_Fortran_COMPILER` options directly.

<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Default</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>CMAKE_C_COMPILER</code></td>
  <td>C Compiler which should be used by CMake</td>
  <td>value of `CC` environment variable at first `cmake` run</td>
</tr>
<tr>
  <td><code>CMAKE_CXX_COMPILER</code></td>
  <td>C++ compiler which should be used by CMake</td>
  <td>
  value of `CXX` environment variable at first `cmake` run
  </td>
</tr>
<tr>
  <td><code>CMAKE_Fortran_COMPILER</code></td>
  <td>C++ compiler which should be used by CMake</td>
  <td>
  value of `FC` environment variable at first `cmake` run
  </td>
</tr>
<tr>
  <td><code>CXX_COMPILER_LAUNCHER</code></td>
  <td>CMake will run this tool and pass the compiler and its arguments to the tool. Some example tools are distcc and ccache.</td>
  <td>
  (empty)
  </td>
</tr>
</tbody>
</table>

### Building with GNU Compilers

```bash
cmake -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ -D CMAKE_Fortran_COMPILER=gfortran ../cmake
```

### Building with Intel Compilers

```bash
cmake -D CMAKE_C_COMPILER=icc -D CMAKE_CXX_COMPILER=icpc -D CMAKE_Fortran_COMPILER=ifort ../cmake
```


### Building with LLVM/Clang Compilers

```bash
cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_Fortran_COMPILER=flang ../cmake
```

## LAMMPS Developer Options


<table>
<thead>
<tr>
  <th>Option</th>
  <th>Description</th>
  <th>Values</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>ENABLE_TESTING</code></td>
  <td>Control wheather to add tests via CTest</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>LAMMPS_TESTING_SOURCE_DIR</code></td>
  <td>Custom location of lammps-testing repository (optional). If not specified it will download it via Git</td>
  <td>
  </td>
</tr>
<tr>
  <td><code>LAMMPS_TESTING_GIT_TAG</code></td>
  <td>If lammps-testing repository is cloned, this is the tag/commit that will be checked out</td>
  <td>
  <dl>
    <dt><code>master</code> (default)</dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_COVERAGE</code></td>
  <td>Enables code coverage support via gcov and adds a gcovr build target to generate a coverage report.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_SANITIZE_ADDRESS</code></td>
  <td>Enables Address Sanitizer support when compiling using GCC or Clang for detecting memory leaks in binaries while running them. See https://clang.llvm.org/docs/AddressSanitizer.html</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_SANITIZE_UNDEFINED</code></td>
  <td>Enables Undefined Behavior Sanitizer support when compiling using GCC or Clang for detecting code that is running into undefined behavior of the language. See https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_SANITIZE_THREAD</code></td>
  <td>Enables Thread Sanitizer support when compiling using GCC or Clang for detecting data races in binaries while running them. See https://clang.llvm.org/docs/ThreadSanitizer.html</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
</tbody>
</table>
