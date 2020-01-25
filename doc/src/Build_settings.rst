Optional build settings
=======================

LAMMPS can be built with several optional settings.  Each sub-section
explain how to do this for building both with CMake and make.

| :ref:`C++11 standard compliance test <cxx11>` when building all of LAMMPS
| :ref:`FFT library <fft>` for use with the :doc:`kspace\_style pppm <kspace_style>` command
| :ref:`Size of LAMMPS data types <size>`
| :ref:`Read or write compressed files <gzip>`
| :ref:`Output of JPG and PNG files <graphics>` via the :doc:`dump image <dump_image>` command
| :ref:`Output of movie files <graphics>` via the :doc:`dump\_movie <dump_image>` command
| :ref:`Memory allocation alignment <align>`
| :ref:`Workaround for long long integers <longlong>`
| :ref:`Error handling exceptions <exceptions>` when using LAMMPS as a library 
| 


----------


.. _cxx11:

C++11 standard compliance test
------------------------------------------

The LAMMPS developers plan to transition to make the C++11 standard the
minimum requirement for compiling LAMMPS.  Currently this only applies to
some packages like KOKKOS while the rest aims to be compatible with the C++98
standard.  Most currently used compilers are compatible with C++11; some need
to set extra flags to switch.  To determine the impact of requiring C++11,
we have added a simple compliance test to the source code, that will cause
the compilation to abort, if C++11 compliance is not available or enabled.
To bypass this check, you need to change a setting in the makefile or
when calling CMake.

**CMake variable**\ :


.. parsed-literal::

   -D DISABLE_CXX11_REQUIREMENT=yes

You can set additional C++ compiler flags (beyond those selected by CMake)
through the CMAKE\_CXX\_FLAGS variable. Example for CentOS 7:


.. parsed-literal::

   -D CMAKE_CXX_FLAGS="-O3 -g -fopenmp -DNDEBUG -std=c++11"

**Makefile.machine setting**\ to bypass the C++11 test and compile in C++98 mode:


.. parsed-literal::

   LMP_INC = -DLAMMPS_CXX98

**Makefile.machine setting**\ to enable the C++11 with older (but not too old) GNU c++ (e.g. on CentOS 7):


.. parsed-literal::

   CCFLAGS = -g -O3 -std=c++11

----------


.. _fft:

FFT library
---------------------

When the KSPACE package is included in a LAMMPS build, the
:doc:`kspace\_style pppm <kspace_style>` command performs 3d FFTs which
require use of an FFT library to compute 1d FFTs.  The KISS FFT
library is included with LAMMPS but other libraries can be faster.
LAMMPS can use them if they are available on your system.

**CMake variables**\ :


.. parsed-literal::

   -D FFT=value              # FFTW3 or MKL or KISS, default is FFTW3 if found, else KISS
   -D FFT_SINGLE=value       # yes or no (default), no = double precision
   -D FFT_PACK=value         # array (default) or pointer or memcpy

.. note::

   The values for the FFT variable must be in upper-case.  This is
   an exception to the rule that all CMake variables can be specified
   with lower-case values.

Usually these settings are all that is needed.  If FFTW3 is selected,
then CMake will try to detect, if threaded FFTW libraries are available
and enable them by default.  This setting is independent of whether
OpenMP threads are enabled and a packages like KOKKOS or USER-OMP is
used.  If CMake cannot detect the FFT library, you can set these variables
to assist:


.. parsed-literal::

   -D FFTW3_INCLUDE_DIRS=path  # path to FFTW3 include files
   -D FFTW3_LIBRARIES=path     # path to FFTW3 libraries
   -D FFT_FFTW_THREADS=on      # enable using threaded FFTW3 libraries
   -D MKL_INCLUDE_DIRS=path    # ditto for Intel MKL library
   -D MKL_LIBRARIES=path

**Makefile.machine settings**\ :


.. parsed-literal::

   FFT_INC = -DFFT_FFTW3         # -DFFT_FFTW3, -DFFT_FFTW (same as -DFFT_FFTW3), -DFFT_MKL, or -DFFT_KISS
                                 # default is KISS if not specified
   FFT_INC = -DFFT_SINGLE        # do not specify for double precision
   FFT_INC = -DFFT_FFTW_THREADS  # enable using threaded FFTW3 libraries
   FFT_INC = -DFFT_PACK_ARRAY    # or -DFFT_PACK_POINTER or -DFFT_PACK_MEMCPY

# default is FFT\_PACK\_ARRAY if not specified


.. parsed-literal::

   FFT_INC =       -I/usr/local/include
   FFT_PATH =      -L/usr/local/lib
   FFT_LIB =       -lfftw3             # FFTW3 double precision
   FFT_LIB =       -lfftw3 -lfftw3_omp # FFTW3 double precision with threads (needs -DFFT_FFTW_THREADS)
   FFT_LIB =       -lfftw3 -lfftw3f    # FFTW3 single precision
   FFT_LIB =       -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  # MKL with Intel compiler
   FFT_LIB =       -lmkl_gf_lp64 -lmkl_sequential -lmkl_core     # MKL with GNU compier

As with CMake, you do not need to set paths in FFT\_INC or FFT\_PATH, if
make can find the FFT header and library files.  You must specify
FFT\_LIB with the appropriate FFT libraries to include in the link.

**CMake and make info**\ :

The `KISS FFT library <http://kissfft.sf.net>`_ is included in the LAMMPS
distribution.  It is portable across all platforms.  Depending on the size
of the FFTs and the number of processors used, the other libraries listed
here can be faster.

However, note that long-range Coulombics are only a portion of the
per-timestep CPU cost, FFTs are only a portion of long-range
Coulombics, and 1d FFTs are only a portion of the FFT cost (parallel
communication can be costly).  A breakdown of these timings is printed
to the screen at the end of a run when using the
:doc:`kspace_style pppm <kspace_style>` command. The :doc:`Run output <Run_output>`
doc page gives more details.  A more detailed (and time consuming)
report of the FFT performance is generated with the
:doc:`kspace_modify fftbench yes <kspace_modify>` command.

FFTW is a fast, portable FFT library that should also work on any
platform and can be faster than the KISS FFT library.  You can
download it from `www.fftw.org <http://www.fftw.org>`_.  LAMMPS requires
version 3.X; the legacy version 2.1.X is no longer supported.

Building FFTW for your box should be as simple as ./configure; make;
make install.  The install command typically requires root privileges
(e.g. invoke it via sudo), unless you specify a local directory with
the "--prefix" option of configure.  Type "./configure --help" to see
various options.

The Intel MKL math library is part of the Intel compiler suite.  It
can be used with the Intel or GNU compiler (see FFT\_LIB setting above).

Performing 3d FFTs in parallel can be time consuming due to data
access and required communication.  This cost can be reduced by
performing single-precision FFTs instead of double precision.  Single
precision means the real and imaginary parts of a complex datum are
4-byte floats.  Double precision means they are 8-byte doubles.  Note
that Fourier transform and related PPPM operations are somewhat less
sensitive to floating point truncation errors and thus the resulting
error is less than the difference in precision. Using the -DFFT\_SINGLE
setting trades off a little accuracy for reduced memory use and
parallel communication costs for transposing 3d FFT data.

When using -DFFT\_SINGLE with FFTW3 you may need to build the FFTW
library a second time with support for single-precision.

For FFTW3, do the following, which should produce the additional
library libfftw3f.a or libfftw3f.so.


.. parsed-literal::

   make clean
   ./configure --enable-single; make; make install

Performing 3d FFTs requires communication to transpose the 3d FFT
grid.  The data packing/unpacking for this can be done in one of 3
modes (ARRAY, POINTER, MEMCPY) as set by the FFT\_PACK syntax above.
Depending on the machine, the size of the FFT grid, the number of
processors used, one option may be slightly faster.  The default is
ARRAY mode.


----------


.. _size:

Size of LAMMPS data types
------------------------------------

LAMMPS has a few integer data types which can be defined as 4-byte or
8-byte integers.  The default setting of "smallbig" is almost always
adequate.

**CMake variable**\ :


.. parsed-literal::

   -D LAMMPS_SIZES=value   # smallbig (default) or bigbig or smallsmall

**Makefile.machine setting**\ :


.. parsed-literal::

   LMP_INC = -DLAMMPS_SMALLBIG    # or -DLAMMPS_BIGBIG or -DLAMMPS_SMALLSMALL

# default is LAMMPS\_SMALLBIG if not specified
**CMake and make info**\ :

The default "smallbig" setting allows for simulations with:

* total atom count = 2\^63 atoms (about 9e18)
* total timesteps = 2\^63 (about 9e18)
* atom IDs = 2\^31 (about 2 billion)
* image flags = roll over at 512

The "bigbig" setting increases the latter two limits.  It allows for:

* total atom count = 2\^63 atoms (about 9e18)
* total timesteps = 2\^63 (about 9e18)
* atom IDs = 2\^63 (about 9e18)
* image flags = roll over at about 1 million (2\^20)

The "smallsmall" setting is only needed if your machine does not
support 8-byte integers.  It allows for:

* total atom count = 2\^31 atoms (about 2 billion)
* total timesteps = 2\^31 (about 2 billion)
* atom IDs = 2\^31 (about 2 billion)
* image flags = roll over at 512 (2\^9)

Atom IDs are not required for atomic systems which do not store bond
topology information, though IDs are enabled by default.  The
:doc:`atom\_modify id no <atom_modify>` command will turn them off.  Atom
IDs are required for molecular systems with bond topology (bonds,
angles, dihedrals, etc).  Thus if you model a molecular system with
more than 2 billion atoms, you need the "bigbig" setting.

Image flags store 3 values per atom which count the number of times an
atom has moved through the periodic box in each dimension.  See the
:doc:`dump <dump>` doc page for a discussion.  If an atom moves through
the periodic box more than this limit, the value will "roll over",
e.g. from 511 to -512, which can cause diagnostics like the
mean-squared displacement, as calculated by the :doc:`compute msd <compute_msd>` command, to be faulty.

Note that the USER-ATC package and the USER-INTEL package are currently
not compatible with the "bigbig" setting. Also, there are limitations
when using the library interface. Some functions with known issues
have been replaced by dummy calls printing a corresponding error rather
than crashing randomly or corrupting data.

Also note that the GPU package requires its lib/gpu library to be
compiled with the same size setting, or the link will fail.  A CMake
build does this automatically.  When building with make, the setting
in whichever lib/gpu/Makefile is used must be the same as above.


----------


.. _graphics:

Output of JPG, PNG, and movie files
--------------------------------------------------

The :doc:`dump image <dump_image>` command has options to output JPEG or
PNG image files.  Likewise the :doc:`dump movie <dump_image>` command
outputs movie files in MPEG format.  Using these options requires the
following settings:

**CMake variables**\ :


.. parsed-literal::

   -D WITH_JPEG=value      # yes or no
                             # default = yes if CMake finds JPEG files, else no
   -D WITH_PNG=value       # yes or no
                             # default = yes if CMake finds PNG and ZLIB files, else no
   -D WITH_FFMPEG=value    # yes or no
                             # default = yes if CMake can find ffmpeg, else no

Usually these settings are all that is needed.  If CMake cannot find
the graphics header, library, executable files, you can set these
variables:


.. parsed-literal::

   -D JPEG_INCLUDE_DIR=path    # path to jpeglib.h header file
   -D JPEG_LIBRARIES=path      # path to libjpeg.a (.so) file
   -D PNG_INCLUDE_DIR=path     # path to png.h header file
   -D PNG_LIBRARIES=path       # path to libpng.a (.so) file
   -D ZLIB_INCLUDE_DIR=path    # path to zlib.h header file
   -D ZLIB_LIBRARIES=path      # path to libz.a (.so) file
   -D FFMPEG_EXECUTABLE=path   # path to ffmpeg executable

**Makefile.machine settings**\ :


.. parsed-literal::

   LMP_INC = -DLAMMPS_JPEG
   LMP_INC = -DLAMMPS_PNG
   LMP_INC = -DLAMMPS_FFMPEG

   JPG_INC = -I/usr/local/include   # path to jpeglib.h, png.h, zlib.h header files if make cannot find them
   JPG_PATH = -L/usr/lib            # paths to libjpeg.a, libpng.a, libz.a (.so) files if make cannot find them
   JPG_LIB = -ljpeg -lpng -lz       # library names

As with CMake, you do not need to set JPG\_INC or JPG\_PATH, if make can
find the graphics header and library files.  You must specify JPG\_LIB
with a list of graphics libraries to include in the link.  You must
insure ffmpeg is in a directory where LAMMPS can find it at runtime,
i.e. a dir in your PATH environment variable.

**CMake and make info**\ :

Using ffmpeg to output movie files requires that your machine
supports the "popen" function in the standard runtime library.

.. note::

   On some clusters with high-speed networks, using the fork()
   library calls (required by popen()) can interfere with the fast
   communication library and lead to simulations using ffmpeg to hang or
   crash.


----------


.. _gzip:

Read or write compressed files
-----------------------------------------

If this option is enabled, large files can be read or written with
gzip compression by several LAMMPS commands, including
:doc:`read\_data <read_data>`, :doc:`rerun <rerun>`, and :doc:`dump <dump>`.

**CMake variables**\ :


.. parsed-literal::

   -D WITH_GZIP=value       # yes or no
                            # default is yes if CMake can find gzip, else no
   -D GZIP_EXECUTABLE=path  # path to gzip executable if CMake cannot find it

**Makefile.machine setting**\ :


.. parsed-literal::

   LMP_INC = -DLAMMPS_GZIP

**CMake and make info**\ :

This option requires that your machine supports the "popen()" function
in the standard runtime library and that a gzip executable can be
found by LAMMPS during a run.

.. note::

   On some clusters with high-speed networks, using the fork()
   library calls (required by popen()) can interfere with the fast
   communication library and lead to simulations using compressed output
   or input to hang or crash. For selected operations, compressed file
   I/O is also available using a compression library instead, which is
   what the :ref:`COMPRESS package <PKG-COMPRESS>` enables.


----------


.. _align:

Memory allocation alignment
---------------------------------------

This setting enables the use of the posix\_memalign() call instead of
malloc() when LAMMPS allocates large chunks or memory.  This can make
vector instructions on CPUs more efficient, if dynamically allocated
memory is aligned on larger-than-default byte boundaries.
On most current systems, the malloc() implementation returns
pointers that are aligned to 16-byte boundaries. Using SSE vector
instructions efficiently, however, requires memory blocks being
aligned on 64-byte boundaries.

**CMake variable**\ :


.. parsed-literal::

   -D LAMMPS_MEMALIGN=value            # 0, 8, 16, 32, 64 (default)

Use a LAMMPS\_MEMALIGN value of 0 to disable using posix\_memalign()
and revert to using the malloc() C-library function instead.  When
compiling LAMMPS for Windows systems, malloc() will always be used
and this setting ignored.

**Makefile.machine setting**\ :


.. parsed-literal::

   LMP_INC = -DLAMMPS_MEMALIGN=value   # 8, 16, 32, 64

Do not set -DLAMMPS\_MEMALIGN, if you want to have memory allocated
with the malloc() function call instead. -DLAMMPS\_MEMALIGN **cannot**
be used on Windows, as it does use different function calls for
allocating aligned memory, that are not compatible with how LAMMPS
manages its dynamical memory.


----------


.. _longlong:

Workaround for long long integers
------------------------------------------------

If your system or MPI version does not recognize "long long" data
types, the following setting will be needed.  It converts "long long"
to a "long" data type, which should be the desired 8-byte integer on
those systems:

**CMake variable**\ :


.. parsed-literal::

   -D LAMMPS_LONGLONG_TO_LONG=value     # yes or no (default)

**Makefile.machine setting**\ :


.. parsed-literal::

   LMP_INC = -DLAMMPS_LONGLONG_TO_LONG


----------


.. _exceptions:

Exception handling when using LAMMPS as a library
------------------------------------------------------------------

This setting is useful when external codes drive LAMMPS as a library.
With this option enabled LAMMPS errors do not kill the caller.
Instead, the call stack is unwound and control returns to the caller,
e.g. to Python.

**CMake variable**\ :


.. parsed-literal::

   -D LAMMPS_EXCEPTIONS=value        # yes or no (default)

**Makefile.machine setting**\ :


.. parsed-literal::

   LMP_INC = -DLAMMPS_EXCEPTIONS


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
