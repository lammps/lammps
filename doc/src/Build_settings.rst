Optional build settings
=======================

LAMMPS can be built with several optional settings.  Each subsection
explains how to do this for building both with CMake and make.

* `C++11 standard compliance`_ when building all of LAMMPS
* `FFT library`_ for use with the :doc:`kspace_style pppm <kspace_style>` command
* `Size of LAMMPS integer types and size limits`_
* `Read or write compressed files`_
* `Output of JPG, PNG, and move files` via the :doc:`dump image <dump_image>` or :doc:`dump movie <dump_image>` commands
* `Memory allocation alignment`_
* `Workaround for long long integers`_
* `Exception handling when using LAMMPS as a library`_ to capture errors
* `Trigger selected floating-point exceptions`_

----------

.. _cxx11:

C++11 standard compliance
------------------------------------------

A C++11 standard compatible compiler is a requirement for compiling LAMMPS.
LAMMPS version 3 March 2020 is the last version compatible with the previous
C++98 standard for the core code and most packages. Most currently used
C++ compilers are compatible with C++11, but some older ones may need extra
flags to enable C++11 compliance.  Example for GNU c++ 4.8.x:

.. code-block:: make

   CCFLAGS = -g -O3 -std=c++11

----------

.. _fft:

FFT library
---------------------

When the KSPACE package is included in a LAMMPS build, the
:doc:`kspace_style pppm <kspace_style>` command performs 3d FFTs which
require use of an FFT library to compute 1d FFTs.  The KISS FFT
library is included with LAMMPS, but other libraries can be faster.
LAMMPS can use them if they are available on your system.

.. versionadded:: TBD

Alternatively, LAMMPS can use the `heFFTe
<https://icl-utk-edu.github.io/heffte/>`_ library for the MPI
communication algorithms, which comes with many optimizations for
special cases, e.g. leveraging available 2D and 3D FFTs in the back end
libraries and better pipelining for packing and communication.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D FFT=value              # FFTW3 or MKL or KISS, default is FFTW3 if found, else KISS
         -D FFT_KOKKOS=value       # FFTW3 or MKL or KISS or CUFFT or HIPFFT, default is KISS
         -D FFT_SINGLE=value       # yes or no (default), no = double precision
         -D FFT_PACK=value         # array (default) or pointer or memcpy
         -D FFT_USE_HEFFTE=value   # yes or no (default), yes links to heFFTe

      .. note::

         When the Kokkos variant of a package is compiled and selected at run time,
         the FFT library selected by the FFT_KOKKOS variable applies. Otherwise,
         the FFT library selected by the FFT variable applies.
         The same FFT settings apply to both. FFT_KOKKOS must be compatible with the
         Kokkos back end - for example, when using the CUDA back end of Kokkos,
         you must use either CUFFT or KISS.

      Usually these settings are all that is needed.  If FFTW3 is
      selected, then CMake will try to detect, if threaded FFTW
      libraries are available and enable them by default.  This setting
      is independent of whether OpenMP threads are enabled and a package
      like KOKKOS or OPENMP is used.  If CMake cannot detect the FFT
      library, you can set these variables to assist:

      .. code-block:: bash

         -D FFTW3_INCLUDE_DIR=path   # path to FFTW3 include files
         -D FFTW3_LIBRARY=path       # path to FFTW3 libraries
         -D FFTW3_OMP_LIBRARY=path   # path to FFTW3 OpenMP wrapper libraries
         -D FFT_FFTW_THREADS=on      # enable using OpenMP threaded FFTW3 libraries
         -D MKL_INCLUDE_DIR=path     # ditto for Intel MKL library
         -D FFT_MKL_THREADS=on       # enable using threaded FFTs with MKL libraries
         -D MKL_LIBRARY=path         # path to MKL libraries
         -D FFT_HEFFTE_BACKEND=value # FFTW or MKL or empty/undefined for the stock heFFTe back end
         -D Heffte_ROOT=path         # path to an existing heFFTe installation

      .. note::

         heFFTe comes with a builtin (= stock) back end for FFTs, i.e. a
         default internal FFT implementation; however, this stock back
         end is intended for testing purposes only and is not optimized
         for production runs.


   .. tab:: Traditional make

      To change the FFT library to be used and its options, you have to edit
      your machine Makefile. Below are examples how the makefile variables
      could be changed.

      .. code-block:: make

         FFT_INC = -DFFT_FFTW3         # -DFFT_FFTW3, -DFFT_FFTW (same as -DFFT_FFTW3), -DFFT_MKL, or -DFFT_KISS
                                       # default is KISS if not specified
         FFT_INC = -DFFT_KOKKOS_CUFFT  # -DFFT_KOKKOS_{FFTW,FFTW3,MKL,CUFFT,HIPFFT,KISS}
                                       # default is KISS if not specified
         FFT_INC = -DFFT_SINGLE        # do not specify for double precision
         FFT_INC = -DFFT_FFTW_THREADS  # enable using threaded FFTW3 libraries
         FFT_INC = -DFFT_MKL_THREADS   # enable using threaded FFTs with MKL libraries
         FFT_INC = -DFFT_PACK_ARRAY    # or -DFFT_PACK_POINTER or -DFFT_PACK_MEMCPY
                                       # default is FFT_PACK_ARRAY if not specified

      .. code-block:: make

         FFT_INC =       -I/usr/local/include
         FFT_PATH =      -L/usr/local/lib
         FFT_LIB =       -lhipfft            # hipFFT either precision
         FFT_LIB =       -lcufft             # cuFFT either precision
         FFT_LIB =       -lfftw3             # FFTW3 double precision
         FFT_LIB =       -lfftw3 -lfftw3_omp # FFTW3 double precision with threads (needs -DFFT_FFTW_THREADS)
         FFT_LIB =       -lfftw3 -lfftw3f    # FFTW3 single precision
         FFT_LIB =       -lmkl_intel_lp64 -lmkl_sequential -lmkl_core   # MKL with Intel compiler, serial interface
         FFT_LIB =       -lmkl_gf_lp64 -lmkl_sequential -lmkl_core      # MKL with GNU compiler, serial interface
         FFT_LIB =       -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core # MKL with Intel compiler, threaded interface
         FFT_LIB =       -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core      # MKL with GNU compiler, threaded interface
         FFT_LIB =       -lmkl_rt            # MKL with automatic runtime selection of interface libs

      As with CMake, you do not need to set paths in ``FFT_INC`` or
      ``FFT_PATH``, if the compiler can find the FFT header and library
      files in its default search path.  You must specify ``FFT_LIB``
      with the appropriate FFT libraries to include in the link.

      Traditional make can also link to heFFTe using an existing installation

      .. code-block:: make

         include <path-to-heffte-installation>/share/heffte/HeffteMakefile.in
         FFT_INC = -DFFT_HEFFTE -DFFT_HEFFTE_FFTW $(heffte_include)
         FFT_PATH =
         FFT_LIB = $(heffte_link) $(heffte_libs)

      The heFFTe install path will contain `HeffteMakefile.in`.
      which will define the `heffte_` include variables needed to link to heFFTe from
      an external project using traditional make.
      The `-DFFT_HEFFTE` is required to switch to using heFFTe, while the optional `-DFFT_HEFFTE_FFTW`
      selects the desired heFFTe back end, e.g., `-DFFT_HEFFTE_FFTW` or `-DFFT_HEFFTE_MKL`,
      omitting the variable will default to the `stock` back end.
      The heFFTe `stock` back end is intended to be used for testing and debugging,
      but is not performance optimized for large scale production runs.

The `KISS FFT library <https://github.com/mborgerding/kissfft>`_ is
included in the LAMMPS distribution.  It is portable across all
platforms.  Depending on the size of the FFTs and the number of
processors used, the other libraries listed here can be faster.

However, note that long-range Coulombics are only a portion of the
per-timestep CPU cost, FFTs are only a portion of long-range Coulombics,
and 1d FFTs are only a portion of the FFT cost (parallel communication
can be costly).  A breakdown of these timings is printed to the screen
at the end of a run when using the :doc:`kspace_style pppm
<kspace_style>` command. The :doc:`Screen and logfile output
<Run_output>` page gives more details.  A more detailed (and time
consuming) report of the FFT performance is generated with the
:doc:`kspace_modify fftbench yes <kspace_modify>` command.

FFTW is a fast, portable FFT library that should also work on any
platform and can be faster than the KISS FFT library.  You can download
it from `www.fftw.org <https://www.fftw.org>`_.  LAMMPS requires version
3.X; the legacy version 2.1.X is no longer supported.

Building FFTW for your box should be as simple as ``./configure; make;
make install``.  The install command typically requires root privileges
(e.g. invoke it via sudo), unless you specify a local directory with
the "--prefix" option of configure.  Type ``./configure --help`` to see
various options.

The Intel MKL math library is part of the Intel compiler suite.  It
can be used with the Intel or GNU compiler (see the ``FFT_LIB`` setting
above).

The cuFFT and hipFFT FFT libraries are packaged with NVIDIA's CUDA and
AMD's HIP installations, respectively. These FFT libraries require the
Kokkos acceleration package to be enabled and the Kokkos back end to be
GPU-resident (i.e., HIP or CUDA).

Performing 3d FFTs in parallel can be time-consuming due to data access
and required communication.  This cost can be reduced by performing
single-precision FFTs instead of double precision.  Single precision
means the real and imaginary parts of a complex datum are 4-byte floats.
Double precision means they are 8-byte doubles.  Note that Fourier
transform and related PPPM operations are somewhat less sensitive to
floating point truncation errors, and thus the resulting error is
generally less than the difference in precision. Using the
``-DFFT_SINGLE`` setting trades off a little accuracy for reduced memory
use and parallel communication costs for transposing 3d FFT data.

When using ``-DFFT_SINGLE`` with FFTW3, you may need to ensure that
the FFTW3 installation includes support for single-precision.

When compiler FFTW3 from source, you can do the following, which should
produce the additional libraries ``libfftw3f.a`` and/or ``libfftw3f.so``\ .

.. code-block:: bash

   make clean
   ./configure --enable-single; make; make install

Performing 3d FFTs requires communication to transpose the 3d FFT
grid.  The data packing/unpacking for this can be done in one of 3
modes (ARRAY, POINTER, MEMCPY) as set by the FFT_PACK syntax above.
Depending on the machine, the size of the FFT grid, the number of
processors used, one option may be slightly faster.  The default is
ARRAY mode.

When using ``-DFFT_HEFFTE`` CMake will first look for an existing
install with hints provided by ``-DHeffte_ROOT``, as recommended by the
CMake standard and note that the name is case sensitive. If CMake cannot
find a heFFTe installation with the correct back end (e.g., FFTW or
MKL), it will attempt to download and build the library automatically.
In this case, LAMMPS CMake will also accept all heFFTe specific
variables listed in the `heFFTe documentation
<https://mkstoyanov.bitbucket.io/heffte/md_doxygen_installation.html>`_
and those variables will be passed into the heFFTe build.

----------

.. _size:

Size of LAMMPS integer types and size limits
--------------------------------------------

LAMMPS uses a few custom integer data types, which can be defined as
either 4-byte (= 32-bit) or 8-byte (= 64-bit) integers at compile time.
This has an impact on the size of a system that can be simulated, or how
large counters can become before "rolling over".  The default setting of
"smallbig" is almost always adequate.

.. tabs::

   .. tab:: CMake build

      With CMake the choice of integer types is made via setting a
      variable during configuration.

      .. code-block:: bash

         -D LAMMPS_SIZES=value   # smallbig (default) or bigbig or smallsmall

      If the variable is not set explicitly, "smallbig" is used.

   .. tab:: Traditional build

      If you want a setting different from the default, you need to edit the
      ``LMP_INC`` variable setting your machine Makefile.

      .. code-block:: make

         LMP_INC = -DLAMMPS_SMALLBIG    # or -DLAMMPS_BIGBIG or -DLAMMPS_SMALLSMALL

      The default setting is ``-DLAMMPS_SMALLBIG`` if nothing is specified

LAMMPS system size restrictions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: auto
   :align: center

   * -
     - smallbig
     - bigbig
     - smallsmall
   * - Total atom count
     - :math:`2^{63}` atoms (= :math:`9.223 \cdot 10^{18}`)
     - :math:`2^{63}` atoms (= :math:`9.223 \cdot 10^{18}`)
     - :math:`2^{31}` atoms (= :math:`2.147 \cdot 10^9`)
   * - Total timesteps
     - :math:`2^{63}` steps (= :math:`9.223 \cdot 10^{18}`)
     - :math:`2^{63}` steps (= :math:`9.223 \cdot 10^{18}`)
     - :math:`2^{31}` steps (= :math:`2.147 \cdot 10^9`)
   * - Atom ID values
     - :math:`1 \le i \le 2^{31} (= 2.147 \cdot 10^9)`
     - :math:`1 \le i \le 2^{63} (= 9.223 \cdot 10^{18})`
     - :math:`1 \le i \le 2^{31} (= 2.147 \cdot 10^9)`
   * - Image flag values
     - :math:`-512 \le i \le 511`
     - :math:`- 1\,048\,576 \le i \le 1\,048\,575`
     - :math:`-512 \le i \le 511`

The "bigbig" setting increases the size of image flags and atom IDs over
"smallbig" and the "smallsmall" setting is only needed if your machine
does not support 64-bit integers or incurs performance penalties when
using them.

These are limits for the core of the LAMMPS code, specific features or
some styles may impose additional limits.  The :ref:`ATC
<PKG-ATC>` package cannot be compiled with the "bigbig" setting.
Also, there are limitations when using the library interface where some
functions with known issues have been replaced by dummy calls printing a
corresponding error message rather than crashing randomly or corrupting
data.

Atom IDs are not required for atomic systems which do not store bond
topology information, though IDs are enabled by default.  The
:doc:`atom_modify id no <atom_modify>` command will turn them off.  Atom
IDs are required for molecular systems with bond topology (bonds,
angles, dihedrals, etc).  Similarly, some force or compute or fix styles
require atom IDs.  Thus, if you model a molecular system or use one of
those styles with more than 2 billion atoms, you need the "bigbig"
setting.

Regardless of the total system size limits, the maximum number of atoms
per MPI rank (local + ghost atoms) is limited to 2 billion for atomic
systems and 500 million for systems with bonds (the additional
restriction is due to using the 2 upper bits of the local atom index
in neighbor lists for storing special bonds info).

Image flags store 3 values per atom in a single integer, which count the
number of times an atom has moved through the periodic box in each
dimension.  See the :doc:`dump <dump>` manual page for a discussion.  If
an atom moves through the periodic box more than this limit, the value
will "roll over", e.g. from 511 to -512, which can cause diagnostics
like the mean-squared displacement, as calculated by the :doc:`compute
msd <compute_msd>` command, to be faulty.

Also note that the GPU package requires its lib/gpu library to be
compiled with the same size setting, or the link will fail.  A CMake
build does this automatically.  When building with make, the setting
in whichever ``lib/gpu/Makefile`` is used must be the same as above.

----------

.. _graphics:

Output of JPG, PNG, and movie files
--------------------------------------------------

The :doc:`dump image <dump_image>` command has options to output JPEG or
PNG image files.  Likewise, the :doc:`dump movie <dump_image>` command
outputs movie files in a variety of movie formats.  Using these options
requires the following settings:

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D WITH_JPEG=value      # yes or no
                                 # default = yes if CMake finds JPEG files, else no
         -D WITH_PNG=value       # yes or no
                                 # default = yes if CMake finds PNG and ZLIB files, else no
         -D WITH_FFMPEG=value    # yes or no
                                 # default = yes if CMake can find ffmpeg, else no

      Usually these settings are all that is needed.  If CMake cannot
      find the graphics header, library, executable files, you can set
      these variables:

      .. code-block:: bash

         -D JPEG_INCLUDE_DIR=path    # path to jpeglib.h header file
         -D JPEG_LIBRARY=path        # path to libjpeg.a (.so) file
         -D PNG_INCLUDE_DIR=path     # path to png.h header file
         -D PNG_LIBRARY=path         # path to libpng.a (.so) file
         -D ZLIB_INCLUDE_DIR=path    # path to zlib.h header file
         -D ZLIB_LIBRARY=path        # path to libz.a (.so) file
         -D FFMPEG_EXECUTABLE=path   # path to ffmpeg executable

   .. tab:: Traditional make

      .. code-block:: make

         LMP_INC = -DLAMMPS_JPEG -DLAMMPS_PNG -DLAMMPS_FFMPEG  <other LMP_INC settings>

         JPG_INC = -I/usr/local/include   # path to jpeglib.h, png.h, zlib.h header files if make cannot find them
         JPG_PATH = -L/usr/lib            # paths to libjpeg.a, libpng.a, libz.a (.so) files if make cannot find them
         JPG_LIB = -ljpeg -lpng -lz       # library names

      As with CMake, you do not need to set ``JPG_INC`` or ``JPG_PATH``,
      if make can find the graphics header and library files in their
      default system locations.  You must specify ``JPG_LIB`` with a
      list of graphics libraries to include in the link.  You must make
      certain that the ffmpeg executable (or ffmpeg.exe on Windows) is
      in a directory where LAMMPS can find it at runtime; that is
      usually a directory list in your ``PATH`` environment variable.

Using ``ffmpeg`` to output movie files requires that your machine
supports the "popen" function in the standard runtime library.

.. note::

   On some clusters with high-speed networks, using the fork()
   library call (required by popen()) can interfere with the fast
   communication library and lead to simulations using ffmpeg to hang or
   crash.

----------

.. _gzip:

Read or write compressed files
-----------------------------------------

If this option is enabled, large files can be read or written with
compression by ``gzip`` or similar tools by several LAMMPS commands,
including :doc:`read_data <read_data>`, :doc:`rerun <rerun>`, and
:doc:`dump <dump>`.  Supported compression tools are currently
``gzip``, ``bzip2``, ``zstd``, and ``lzma``.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D WITH_GZIP=value       # yes or no
                                  # default is yes if CMake can find the gzip program, else no

   .. tab:: Traditional make

      .. code-block:: make

         LMP_INC = -DLAMMPS_GZIP   <other LMP_INC settings>

This option requires that your operating system fully supports the
"popen()" function in the standard runtime library and that a ``gzip``
or other executable can be found by LAMMPS in the standard search path
during a run.

.. note::

   On clusters with high-speed networks, using the "fork()" library call
   (required by "popen()") can interfere with the fast communication
   library and lead to simulations using compressed output or input to
   hang or crash. For selected operations, compressed file I/O is also
   available using a compression library instead, which is what the
   :ref:`COMPRESS package <PKG-COMPRESS>` enables.

----------

.. _align:

Memory allocation alignment
---------------------------------------

This setting enables the use of the "posix_memalign()" call instead of
"malloc()" when LAMMPS allocates large chunks of memory.  Vector
instructions on CPUs may become more efficient, if dynamically allocated
memory is aligned on larger-than-default byte boundaries.  On most
current operating systems, the "malloc()" implementation returns
pointers that are aligned to 16-byte boundaries. Using SSE vector
instructions efficiently, however, requires memory blocks being aligned
on 64-byte boundaries.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D LAMMPS_MEMALIGN=value            # 0, 8, 16, 32, 64 (default)

      Use a ``LAMMPS_MEMALIGN`` value of 0 to disable using
      "posix_memalign()" and revert to using the "malloc()" C-library
      function instead.  When compiling LAMMPS for Windows systems,
      "malloc()" will always be used and this setting is ignored.

   .. tab:: Traditional make

      .. code-block:: make

         LMP_INC = -DLAMMPS_MEMALIGN=value   # 8, 16, 32, 64

      Do not set ``-DLAMMPS_MEMALIGN``, if you want to have memory
      allocated with the "malloc()" function call
      instead. ``-DLAMMPS_MEMALIGN`` **cannot** be used on Windows, as
      Windows different function calls with different semantics for
      allocating aligned memory, that are not compatible with how LAMMPS
      manages its dynamical memory.

----------

.. _longlong:

Workaround for long long integers
---------------------------------

If your system or MPI version does not recognize "long long" data
types, the following setting will be needed.  It converts "long long"
to a "long" data type, which should be the desired 8-byte integer on
those systems:

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D LAMMPS_LONGLONG_TO_LONG=value     # yes or no (default)

   .. tab:: Traditional make

      .. code-block:: make

         LMP_INC = -DLAMMPS_LONGLONG_TO_LONG  <other LMP_INC settings>

----------

.. _exceptions:

Exception handling when using LAMMPS as a library
-------------------------------------------------

LAMMPS errors do not kill the calling code, but throw an exception.  In
the C-library interface, the call stack is unwound and control returns
to the caller, e.g. to Python or a code that is coupled to LAMMPS. The
error status can then be queried.  When using C++ directly, the calling
code has to be set up to *catch* exceptions thrown from within LAMMPS.

.. note::

   When LAMMPS is running in parallel, it is not always possible to
   cleanly recover from an exception since not all parallel ranks may
   throw an exception and thus other MPI ranks may get stuck waiting for
   messages from the ones with errors.

----------

.. _trap_fpe:

Trigger selected floating-point exceptions
------------------------------------------

Many kinds of CPUs have the capability to detect when a calculation
results in an invalid math operation, like a division by zero or calling
the square root with a negative argument.  The default behavior on
most operating systems is to continue and have values for ``NaN`` (= not
a number) or ``Inf`` (= infinity).  This allows software to detect and
recover from such conditions.  This behavior can be changed, however,
often through use of compiler flags.  On Linux systems (or more general
on systems using the GNU C library), these so-called floating-point traps
can also be selectively enabled through library calls.  LAMMPS supports
that by setting the ``-DLAMMPS_TRAP_FPE`` pre-processor define.  As it is
done in the ``main()`` function, this applies only to the standalone
executable, not the library.

.. tabs::

   .. tab:: CMake build

      .. code-block:: bash

         -D CMAKE_TUNE_FLAGS=-DLAMMPS_TRAP_FPE

   .. tab:: Traditional make

      .. code-block:: make

         LMP_INC = -DLAMMPS_TRAP_FPE  <other LMP_INC settings>

After compilation with this flag set, the LAMMPS executable will stop
and produce a core dump when a division by zero, overflow, illegal math
function argument or other invalid floating point operation is encountered.
