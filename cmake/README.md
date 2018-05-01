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

## Quick Start for the Impatient
If you want to skip ahead and just run the compilation using `cmake`, please
find a minimal example below. Together with the options reference below, this
should get you started.

```bash
git clone https://github.com/lammps/lammps.git
mkdir lammps/build
cd lammps/build
cmake ../cmake [-DOPTION_A=VALUE_A -DOPTION_B=VALUE_B ...]
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

# Reference

## Common CMAKE Configuration Options


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
  </td>
</tr>
<tr>
  <td><code>CMAKE_BUILD_TYPE</code></td>
  <td>Controls if debugging symbols are added to the generated binaries</td>
  <td>
  <dl>
  <dt><code>Release</code> (default)</dt>
  <dt><code>Debug</code></dt>
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
  <td>allows appending a machine suffix to the generate LAMMPS binary</td>
  <td>
  <dl>
    <dt>*none*  (default)</dt>
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
  <td><code>LAMMPS_LONGLONG_TO_LONG</code></td>
  <td>Workaround if your system or MPI version does not recognize <code>long long</code> data types</td>
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
  <td><code>ENABLE_MPI</code></td>
  <td>control whether to build LAMMPS with MPI support. This will look for
  `mpicxx` in your path and use this MPI implementation.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_OPT</code></td>
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
  <td><code>ENABLE_USER-OMP</code></td>
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
  <td><code>ENABLE_USER-INTEL</code></td>
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
  <td><code>ENABLE_GPU</code></td>
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
  <td><code>ENABLE_KOKKOS</code></td>
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
  <td><code>ENABLE_ALL</code></td>
  <td>Enable all default packages</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_ASPHERE</code></td>
  <td>Computes, time-integration fixes, and pair styles for aspherical particle models including ellipsoids, 2d lines, and 3d triangles.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_BODY</code></td>
  <td>Body-style particles with internal structure. Computes, time-integration fixes, pair styles, as well as the body styles themselves.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_CLASS2</code></td>
  <td>Bond, angle, dihedral, improper, and pair styles for the COMPASS CLASS2 molecular force field.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_COLLOID</code></td>
  <td>Coarse-grained finite-size colloidal particles. Pair styles and fix wall styles for colloidal interactions. Includes the Fast Lubrication Dynamics (FLD) method for hydrodynamic interactions, which is a simplified approximation to Stokesian dynamics.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_COMPRESS</code></td>
  <td>Compressed output of dump files via the zlib compression library, using dump styles with a “gz” in their style name.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_CORESHELL</code></td>
  <td>Compute and pair styles that implement the adiabatic core/shell model for polarizability. The pair styles augment Born, Buckingham, and Lennard-Jones styles with core/shell capabilities. The compute temp/cs command calculates the temperature of a system with core/shell particles.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_DIPOLE</code></td>
  <td>An atom style and several pair styles for point dipole models with short-range or long-range interactions.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_GRANULAR</code></td>
  <td>Pair styles and fixes for finite-size granular particles, which interact with each other and boundaries via frictional and dissipative potentials.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_KSPACE</code></td>
  <td>A variety of long-range Coulombic solvers, as well as pair styles which compute the corresponding short-range pairwise Coulombic interactions. These include Ewald, particle-particle particle-mesh (PPPM), and multilevel summation method (MSM) solvers.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_MANYBODY</code></td>
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
  <td><code>ENABLE_MC</code></td>
  <td>
  Several fixes and a pair style that have Monte Carlo (MC) or MC-like
  attributes. These include fixes for creating, breaking, and swapping bonds,
  for performing atomic swaps, and performing grand-canonical MC (GCMC) in
  conjuction with dynamics.
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_MEAM</code></td>
  <td>
<p>A pair style for the modified embedded atom (MEAM) potential.</p>

<p><strong>Please note that the MEAM package has been superseded by the USER-MEAMC package,
which is a direct translation of the MEAM package to C++. USER-MEAMC contains
additional optimizations making it run faster than MEAM on most machines, while
providing the identical features and USER interface.</strong></p>
  </td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_MISC</code></td>
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
  <td><code>ENABLE_MOLECULE</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_PERI</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_QEQ</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_REAX</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_REPLICA</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_RIGID</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_SHOCK</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_SNAP</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_SRD</code></td>
  <td></td>
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
  <td><code>ENABLE_KIM</code></td>
  <td>A <code>pair_style kim</code> command which is a wrapper on the Knowledge Base for Interatomic Models (KIM) repository of interatomic potentials, enabling any of them to be used in LAMMPS simulations.</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_PYTHON</code></td>
  <td>Enable all default packages</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_MSCG</code></td>
  <td>Enable all default packages</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_MPIIO</code></td>
  <td>Enable all default packages</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_POEMS</code></td>
  <td>Enable all default packages</td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_LATTE</code></td>
  <td>Enable all default packages</td>
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
  <td><code>ENABLE_USER-ATC</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-AWPMD</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-CGDNA</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-MESO</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-CGSDK</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-COLVARS</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-DIFFRACTION</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-DPD</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-DRUDE</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-EFF</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-FEP</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-H5MD</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-LB</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-MANIFOLD</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-MEAMC</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-MGPT</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-MISC</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-MOFFF</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-MOLFILE</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-NETCDF</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-PHONON</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-QTB</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-REAXC</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-SMD</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-SMTBQ</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-SPH</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-UEF</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-VTK</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-QUIP</code></td>
  <td></td>
  <td>
  <dl>
    <dt><code>off</code> (default)</dt>
    <dt><code>on</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>ENABLE_USER-QMMM</code></td>
  <td></td>
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
    <p>If either MKL or FFTW is selected <code>cmake</code> will try to locate these libraries automatically. To control which one should be used please see the options below for each FFT library.</p>
  </td>
  <td>
  <dl>
    <dt><code>KISSFFT</code></dt>
    <dt><code>FFTW3</code></dt>
    <dt><code>FFTW2</code></dt>
    <dt><code>MKL</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>PACK_ARRAY</code></td>
  <td>Optimization for FFT</td>
  <td>
  <dl>
    <dt><code>PACK_ARRAY</code></dt>
    <dt><code>PACK_POINTER</code></dt>
    <dt><code>PACK_MEMCPY</code></dt>
  </dl>
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
  <td>Optimization for FFT</td>
  <td>
  </td>
</tr>
<tr>
  <td><code>MKL_LIBRARIES</code></td>
  <td>Optimization for FFT</td>
  <td>
  </td>
</tr>
</tbody>
</table>

TODO static vs dynamic linking

### FFTW2

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
  <td><code>FFTW2_INCLUDE_DIRS</code></td>
  <td>Optimization for FFT</td>
  <td>
  </td>
</tr>
<tr>
  <td><code>FFTW2_LIBRARIES</code></td>
  <td>Optimization for FFT</td>
  <td>
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
  <td>Optimization for FFT</td>
  <td>
  </td>
</tr>
<tr>
  <td><code>FFTW3_LIBRARIES</code></td>
  <td>Optimization for FFT</td>
  <td>
  </td>
</tr>
</tbody>
</table>

### LAPACK
TODO

### PYTHON Package


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
    <dt><code>OpenCL</code> (default)</dt>
    <dt><code>CUDA</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>GPU_PREC</code></td>
  <td>Precision size used by GPU package kernels</td>
  <td>
  <dl>
    <dt><code>SINGLE_DOUBLE</code></dt>
    <dt><code>SINGLE_SINGLE</code></dt>
    <dt><code>DOUBLE_DOUBLE</code></dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>OCL_TUNE</code> (OpenCL only)</td>
  <td>Tuning target for OpenCL driver code</td>
  <td>
  <dl>
    <dt><code>GENERIC</code> (default)</dt>
    <dt><code>INTEL</code> (Intel CPU)</dt>
    <dt><code>PHI</code> (Intel Xeon Phi)</dt>
    <dt><code>FERMI</code> (NVIDIA)</dt>
    <dt><code>KEPLER</code> (NVIDIA)</dt>
    <dt><code>CYPRESS</code> (AMD)</dt>
  </dl>
  </td>
</tr>
<tr>
  <td><code>GPU_ARCH</code> (CUDA only)</td>
  <td>CUDA SM architecture targeted by GPU package</td>
  <td>
  <dl>
    <dt><code>sm20</code> (Fermi)</dt>
    <dt><code>sm30</code> (Kepler)</dt>
    <dt><code>sm50</code> (Maxwell)</dt>
    <dt><code>sm60</code> (Pascal)</dt>
    <dt><code>sm70</code> (Volta)</dt>
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
</tbody>
</table>

### VORONOI Package

TODO

### USER-SMD Package

Requires a Eigen3 installation

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
  <td><code>EIGEN3_INCLUDE_DIR</code></td>
  <td></td>
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
  <td><code>GZIP_EXECUTABLE</code></td>
  <td></td>
  <td>
  </td>
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
  <td><code>FFMPEG_EXECUTABLE</code></td>
  <td></td>
  <td>
  </td>
</tr>
</tbody>
</table>


## Compilers

By default, `cmake` will use your environment C/C++/Fortran compilers for a build. It uses the `CC`, `CXX` and `FC` environment variables to detect which compilers should be used. However, these values
will be cached after the first run of `cmake`. Subsequent runs of `cmake` will ignore changes in these environment variables. To ensure the correct values are used you avoid the cache by setting the `CMAKE_C_COMPILER`, `CMAKE_CXX_COMPILER`, `CMAKE_Fortran_COMPILER` options directly.

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
</tbody>
</table>

### Building with GNU Compilers

```bash
cmake ../cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran
```

### Building with Intel Compilers

```bash
cmake ../cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort
```


### Building with LLVM/Clang Compilers

```bash
cmake ../cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_Fortran_COMPILER=flang
```


## Examples
