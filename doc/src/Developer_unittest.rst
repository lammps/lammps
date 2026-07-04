Adding tests for unit testing
-----------------------------

.. contents::
   :local:

------------

This section discusses adding or expanding tests for the unit test
infrastructure included into the LAMMPS source code distribution.
Unlike example inputs, unit tests focus on testing the "local" behavior
of individual features, tend to run fast, and should be set up to cover
as much of the added code as possible.  When contributing code to the
distribution, the LAMMPS developers will appreciate if additions to the
integrated unit test facility are included.

Given the complex nature of MD simulations where many operations can
only be performed when suitable "real" simulation environment has been
set up, not all tests will be unit tests in the strict definition of
the term.  They are rather executed on a more abstract level by issuing
LAMMPS script commands and then inspecting the changes to the internal
data.  For some classes of tests, generic test programs have been
written that can be applied to parts of LAMMPS that use the same
interface (via polymorphism) and those are driven by input files, so
tests can be added by simply adding more of those input files.  Those
tests should be seen more as a hybrid between unit and regression tests.

When adding tests it is recommended to also :ref:`enable support for
code coverage reporting <testing>`, and study the coverage reports
so that it is possible to monitor which parts of the code of a given
file are executed during the tests and which tests would need to be
added to increase the coverage.

The tests are grouped into categories and corresponding folders.
The following sections describe how the tests are implemented and
executed in those categories with increasing complexity of tests
and implementation.


Tests for utility functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

These tests are driven by programs in the ``unittest/utils`` folder
and most closely resemble conventional unit tests. There is one test
program for each namespace or group of classes or file. The naming
convention for the sources and executables is that they start with
with ``test_``.  The following sources and groups of tests are currently
available:

.. list-table::
   :header-rows: 1
   :widths: 32 18 50
   :align: left

   * - File name:
     - Test name:
     - Description:
   * - ``test_argutils.cpp``
     - ArgInfo
     - Tests for ``ArgInfo`` class used by LAMMPS
   * - ``test_fmtlib.cpp``
     - FmtLib
     - Tests for ``{fmt}`` or ``std::format`` functions used by LAMMPS
   * - ``test_math_eigen_impl.cpp``
     - MathEigen
     - Tests for ``MathEigen::`` classes and functions
   * - ``test_mempool.cpp``
     - MemPool
     - Tests for :cpp:class:`MyPage <LAMMPS_NS::MyPage>` and :cpp:class:`MyPoolChunk <LAMMPS_NS::MyPoolChunk>`
   * - ``test_tokenizer.cpp``
     - Tokenizer
     - Tests for :cpp:class:`Tokenizer <LAMMPS_NS::Tokenizer>` and :cpp:class:`ValueTokenizer <LAMMPS_NS::ValueTokenizer>`
   * - ``test_utils.cpp``
     - Utils
     - Tests for ``utils::`` :doc:`functions <Developer_utils>`
   * - ``test_fft3d.cpp``
     - FFT3D
     - Tests for standard FFT3d wrapper (KISS, FFTW3, MKL, NVPL)
   * - ``test_fft3d_kokkos.cpp``
     - FFT3DKokkos
     - Tests for KOKKOS FFT3d wrapper (CPU and GPU back ends)


To add tests either an existing source file needs to be modified or a
new source file needs to be added to the distribution and enabled for
testing.  To add a new file suitable CMake script code needs to be added
to the ``CMakeLists.txt`` file in the ``unittest/utils`` folder.  Example:

.. code-block:: cmake

   add_executable(test_tokenizer test_tokenizer.cpp)
   target_link_libraries(test_tokenizer PRIVATE lammps GTest::GMockMain GTest::GMock GTest::GTest)
   add_test(Tokenizer test_tokenizer)

This adds instructions to build the ``test_tokenizer`` executable from
``test_tokenizer.cpp`` and links it with the GoogleTest libraries and the
LAMMPS library as well as it uses the ``main()`` function from the
GoogleMock library of GoogleTest.  The third line registers the executable
as a test program to be run from ``ctest`` under the name ``Tokenizer``.

The test executable itself will execute multiple individual tests
through the GoogleTest framework. In this case each test consists of
creating a tokenizer class instance with a given string and explicit or
default separator choice, and then executing member functions of the
class and comparing their results with expected values. A few examples:

.. code-block:: c++

   TEST(Tokenizer, empty_string)
   {
       Tokenizer t("", " ");
       ASSERT_EQ(t.count(), 0);
   }

   TEST(Tokenizer, two_words)
   {
       Tokenizer t("test word", " ");
       ASSERT_EQ(t.count(), 2);
   }

   TEST(Tokenizer, default_separators)
   {
       Tokenizer t(" \r\n test \t word \f");
       ASSERT_THAT(t.next(), Eq("test"));
       ASSERT_THAT(t.next(), Eq("word"));
       ASSERT_EQ(t.count(), 2);
   }

Each of these TEST functions will become an individual
test run by the test program. When using the ``ctest``
command as a front end to run the tests, their output
will be suppressed and only a summary printed, but adding
the '-V' option will then produce output from the tests
above like the following:

.. code-block:: console

   [...]
   1: [ RUN      ] Tokenizer.empty_string
   1: [       OK ] Tokenizer.empty_string (0 ms)
   1: [ RUN      ] Tokenizer.two_words
   1: [       OK ] Tokenizer.two_words (0 ms)
   1: [ RUN      ] Tokenizer.default_separators
   1: [       OK ] Tokenizer.default_separators (0 ms)
   [...]

The MathEigen test collection has been adapted from a standalone test
and does not use the GoogleTest framework and thus not representative.
The other test sources, however, can serve as guiding examples for
additional tests.

FFT Testing Infrastructure
""""""""""""""""""""""""""

.. versionadded:: 10Dec2025

The FFT tests (``test_fft3d.cpp`` and ``test_fft3d_kokkos.cpp``)
validate the LAMMPS FFT wrapper implementations for both standard (CPU)
and KOKKOS (CPU/GPU) back ends.  These tests require the KSPACE package
and use specialized helper utilities to ensure FFT correctness across
different library back ends (KISS FFT, FFTW3, MKL, NVPL, cuFFT, hipFFT,
etc.).

**Building and Running FFT Tests:**

The FFT tests are automatically enabled when ``ENABLE_TESTING=ON`` and
``PKG_KSPACE=ON`` are set during CMake configuration. For KOKKOS FFT tests,
``PKG_KOKKOS=ON`` is also required.

Run only FFT tests using the ``ctest`` command of the CMake software:

.. code-block:: bash

   ctest -R FFT3D          # Run all tests with FFT3D in their name
   ctest -R FFT3D -V       # Same as above but with verbose output
   ctest -L fft            # Run all tests labeled with 'fft'

Tests automatically skip configurations requiring libraries or back ends
not available in the current build (e.g., FFTW3, MPI, CUDA).

**FFT Test Helper Header:**

The testing infrastructure uses ``fft_test_helpers.h`` which contains
test data generators, validators, and utilities.

For runtime configuration detection, tests use the existing ``Info``
class API (``Info::has_package()``, ``Info::has_accelerator_feature()``,
etc.).

The ``fft_test_helpers.h`` header provides three main namespaces:

1. **FFTTestHelpers** - utility functions:
   ``FFTBuffer`` (RAII wrapper), ``idx3d()`` (index conversion),
   ``scaled_tolerance()`` (grid-size-dependent precision)

2. **FFTTestData** - test data generators:
   ``DeltaFunctionGenerator``, ``ConstantGenerator``, ``SineWaveGenerator``,
   ``GaussianGenerator``, ``RandomComplexGenerator``, ``MixedModesGenerator``

3. **FFTValidation** - validators:
   ``RoundTripValidator``, ``KnownAnswerValidator``, ``ParsevalValidator``,
   ``HermitianSymmetryValidator``, ``LinearityValidator``

**Example Test:**

.. code-block:: c++

   TEST_F(FFT3DTest, RoundTrip_32x32x32) {
       FFTBuffer original(32, 32, 32), fft_result(32, 32, 32), recovered(32, 32, 32);

       GaussianGenerator generator(2.0);
       generator.generate(original.data(), 32, 32, 32);

       fft->compute(original.data(), fft_result.data(), FFT3d::FORWARD);
       fft->compute(fft_result.data(), recovered.data(), FFT3d::BACKWARD);

       RoundTripValidator validator(original.data(), recovered.data(), 32, 32, 32,
                                     scaled_tolerance(ROUNDTRIP_TOLERANCE, 32, 32, 32));
       EXPECT_TRUE(validator.validate());
   }

**Precision and Tolerances:**

FFT tests use precision-aware tolerances that automatically adjust based
on floating-point precision (``-D FFT_SINGLE=ON`` vs ``-D
FFT_SINGLE=off``), grid size, and accelerator back end.  Base tolerances
(``ROUNDTRIP_TOLERANCE``, ``PARSEVAL_TOLERANCE``, etc.)  are defined in
``fft_test_helpers.h``.  Use ``scaled_tolerance()`` to adjust for grid
size effects.

Tests for individual LAMMPS commands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tests ``unittest/commands`` are a bit more complex as they require
to first create a :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class instance
and then use the :doc:`C++ API <Cplusplus>` to pass individual commands
to that LAMMPS instance.  For that reason these tests use a GoogleTest
"test fixture", i.e. a class derived from ``testing::Test`` that will
create (and delete) the required LAMMPS class instance for each set of
tests in a ``TEST_F()`` function.  Please see the individual source files
for different examples of setting up suitable test fixtures.  Here is an
example for implementing a test using a fixture by first checking the
default value and then issuing LAMMPS commands and checking whether they
have the desired effect:

.. code-block:: c++

   TEST_F(SimpleCommandsTest, ResetTimestep)
   {
       ASSERT_EQ(lmp->update->ntimestep, 0);

       BEGIN_HIDE_OUTPUT();
       command("reset_timestep 10");
       END_HIDE_OUTPUT();
       ASSERT_EQ(lmp->update->ntimestep, 10);

       BEGIN_HIDE_OUTPUT();
       command("reset_timestep 0");
       END_HIDE_OUTPUT();
       ASSERT_EQ(lmp->update->ntimestep, 0);

       TEST_FAILURE(".*ERROR: Timestep must be >= 0.*", command("reset_timestep -10"););
       TEST_FAILURE(".*ERROR: Illegal reset_timestep .*", command("reset_timestep"););
       TEST_FAILURE(".*ERROR: Illegal reset_timestep .*", command("reset_timestep 10 10"););
       TEST_FAILURE(".*ERROR: Expected integer .*", command("reset_timestep xxx"););
   }

Please note the use of the ``BEGIN_HIDE_OUTPUT`` and ``END_HIDE_OUTPUT``
functions that will capture output from running LAMMPS.  This is normally
discarded but by setting the verbose flag (via setting the ``TEST_ARGS``
environment variable, ``TEST_ARGS=-v``) it can be printed and used to
understand why tests fail unexpectedly.

The specifics of so-called "death tests", i.e. conditions where LAMMPS
should fail and throw an exception, are implemented in the
``TEST_FAILURE()`` macro. These tests operate by capturing the screen
output when executing the failing command and then comparing that with a
provided regular expression string pattern.  Example:

.. code-block:: c++

   TEST_F(SimpleCommandsTest, UnknownCommand)
   {
       TEST_FAILURE(".*ERROR: Unknown command.*", lmp->input->one("XXX one two"););
   }

The following test programs are currently available:

.. list-table::
   :header-rows: 1
   :widths: auto
   :align: left

   * - File name:
     - Test name:
     - Description:
   * - ``test_simple_commands.cpp``
     - SimpleCommands
     - Tests for LAMMPS commands that do not require a box
   * - ``test_lattice_region.cpp``
     - LatticeRegion
     - Tests to validate the :doc:`lattice <lattice>` and :doc:`region <region>` commands
   * - ``test_groups.cpp``
     - GroupTest
     - Tests to validate the :doc:`group <group>` command
   * - ``test_variables.cpp``
     - VariableTest
     - Tests to validate the :doc:`variable <variable>` command
   * - ``test_kim_commands.cpp``
     - KimCommands
     - Tests for several commands from the :ref:`KIM package <PKG-KIM>`
   * - ``test_reset_atoms.cpp``
     - ResetAtoms
     - Tests to validate the :doc:`reset_atoms <reset_atoms>` sub-commands


Tests for the C-style library interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tests for validating the LAMMPS C-style library interface are in the
``unittest/c-library`` folder.  They text either utility functions or
LAMMPS commands, but use the functions implemented in
``src/library.cpp`` as much as possible.  There may be some overlap with
other tests as far as the LAMMPS functionality is concerned, but the
focus is on testing the C-style library API.  The tests are distributed
over multiple test programs which try to match the grouping of the
functions in the source code and :ref:`in the manual <lammps_c_api>`.

This group of tests also includes tests invoking LAMMPS in parallel
through the library interface, provided that LAMMPS was compiled with
MPI support.  These include tests where LAMMPS is run in multi-partition
mode or only on a subset of the MPI world communicator.  The CMake
script code for adding this kind of test looks like this:

.. code-block:: cmake

   if (BUILD_MPI)
     add_executable(test_library_mpi test_library_mpi.cpp)
     target_link_libraries(test_library_mpi PRIVATE lammps GTest::GTest GTest::GMock)
     target_compile_definitions(test_library_mpi PRIVATE ${TEST_CONFIG_DEFS})
     add_mpi_test(NAME LibraryMPI NUM_PROCS 4 COMMAND $<TARGET_FILE:test_library_mpi>)
   endif()

Note the custom function ``add_mpi_test()`` which adapts how ``ctest``
will execute the test so it is launched in parallel (with 4 MPI ranks).

Tests for the Python module and package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``unittest/python`` folder contains primarily tests for classes and
functions in the LAMMPS python module but also for commands in the
PYTHON package.  These tests are only enabled, if the necessary
prerequisites are detected or enabled during configuration and
compilation of LAMMPS (shared library build enabled, Python interpreter
found, Python development files found).

The Python tests are implemented using the ``unittest`` standard Python
module and split into multiple files with similar categories as the
tests for the C-style library interface.

Tests for the Fortran interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tests for using the Fortran module are in the ``unittest/fortran``
folder.  Since they are also using the GoogleTest library, they require
test wrappers written in C++ that will call fortran functions with a C
function interface through ISO_C_BINDINGS which will in turn call the
functions in the LAMMPS Fortran module.

Tests for the C++-style library interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tests in the ``unittest/cplusplus`` folder are somewhat similar to
the tests for the C-style library interface, but do not need to test the
convenience and utility functions that are only available through the
C-style library interface.  Instead they focus on the more generic
features that are used in LAMMPS internally.  This part of the unit
tests is currently still mostly in the planning stage.

Tests for reading and writing file formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``unittest/formats`` folder contains test programs for reading and
writing files like data files, restart files, potential files or dump
files.  This covers simple things like the file i/o convenience
functions in the ``utils::`` namespace to complex tests of atom styles
where creating and deleting of atoms with different properties is tested
in different ways and through script commands or reading and writing of
data or restart files.

Tests for styles computing or modifying forces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are tests common configurations for pair styles, bond styles,
angle styles, kspace styles and certain fix styles.  Those are tests
driven by some test executables build from sources in the
``unittest/force-styles`` folder and use LAMMPS input template and data
files as well as input files in YAML format from the
``unittest/force-styles/tests`` folder. The YAML file names have to
follow some naming conventions so they get associated with the test
programs and categorized and listed with canonical names in the list
of tests as displayed by ``ctest -N``.  If you add a new YAML file,
you need to re-run CMake to update the corresponding list of tests.

A minimal YAML file for a (molecular) pair style test will looks
something like the following (see ``mol-pair-zero.yaml``):

.. code-block:: yaml

   ---
   lammps_version: 24 Aug 2020
   date_generated: Tue Sep 15 09:44:21 202
   epsilon: 1e-14
   prerequisites: ! |
     atom full
     pair zero
   pre_commands: ! ""
   post_commands: ! ""
   input_file: in.fourmol
   pair_style: zero 8.0
   pair_coeff: ! |
     * *
   extract: ! ""
   natoms: 29
   init_vdwl: 0
   init_coul: 0

   [...]

The following table describes the available keys and their purpose.  The
``Tester`` column lists which test program(s) use each key.  ``all`` means
every force-style tester (``test_pair_style``, ``test_bond_style``,
``test_angle_style``, ``test_dihedral_style``, ``test_improper_style``, and
``test_fix_timestep``), and ``bonded interaction tests`` means the
``test_bond_style``, ``test_angle_style``, ``test_dihedral_style``, and
``test_improper_style`` testers.  The ``Required`` column indicates whether
a valid YAML file (for the listed tester) must contain the key (``yes``) or
whether it may be omitted (``no``).  The reference generator always writes
the required keys; the optional keys are either metadata added by hand
(``skip_tests``, ``tags``) or reference data that is only emitted when the
tested style actually provides it (e.g. ``global_scalar`` only for a fix
with scalar output):

.. list-table::
   :header-rows: 1

   * - Key:
     - Tester:
     - Required:
     - Description:
   * - lammps_version
     - all
     - yes
     - LAMMPS version used to last update the reference data
   * - date_generated
     - all
     - yes
     - date when the file was last updated
   * - epsilon
     - all
     - yes
     - base value for the relative precision required for tests to pass
   * - skip_tests
     - all
     - no
     - request to skip the indicated test fixtures (see table below)
   * - tags
     - all
     - no
     - used to classify tests and to adjust behavior of test fixtures (see table below)
   * - prerequisites
     - all
     - yes
     - list of style kind / style name pairs required to run the test
   * - pre_commands
     - all
     - yes
     - LAMMPS commands to be executed before the input template file is read
   * - post_commands
     - all
     - yes
     - LAMMPS commands to be executed right before the actual tests
   * - input_file
     - all
     - yes
     - LAMMPS input file template
   * - input_coeffs
     - test_fix_timestep
     - yes
     - file with the force-field and group setup commands applied after the input template
   * - natoms
     - all
     - yes
     - number of atoms in the input file template
   * - pair_style
     - test_pair_style
     - yes
     - arguments to the pair_style command to be tested
   * - pair_coeff
     - test_pair_style
     - yes
     - list of pair_coeff arguments to set parameters for the input template
   * - init_vdwl
     - test_pair_style
     - yes
     - non-Coulomb pair energy after "run 0"
   * - init_coul
     - test_pair_style
     - yes
     - Coulomb pair energy after "run 0"
   * - run_vdwl
     - test_pair_style
     - yes
     - non-Coulomb pair energy after "run 4"
   * - run_coul
     - test_pair_style
     - yes
     - Coulomb pair energy after "run 4"
   * - bond_style
     - test_bond_style
     - yes
     - arguments to the bond_style command to be tested
   * - bond_coeff
     - test_bond_style
     - yes
     - list of bond_coeff arguments to set parameters
   * - angle_style
     - test_angle_style
     - yes
     - arguments to the angle_style command to be tested
   * - angle_coeff
     - test_angle_style
     - yes
     - list of angle_coeff arguments to set parameters
   * - dihedral_style
     - test_dihedral_style
     - yes
     - arguments to the dihedral_style command to be tested
   * - dihedral_coeff
     - test_dihedral_style
     - yes
     - list of dihedral_coeff arguments to set parameters
   * - improper_style
     - test_improper_style
     - yes
     - arguments to the improper_style command to be tested
   * - improper_coeff
     - test_improper_style
     - yes
     - list of improper_coeff arguments to set parameters
   * - init_energy
     - bonded interaction tests
     - yes
     - bonded interaction energy after "run 0"
   * - run_energy
     - bonded interaction tests
     - yes
     - bonded interaction energy after "run 4"
   * - equilibrium
     - test_bond_style
     - yes
     - equilibrium distance for each type
   * - equilibrium
     - test_angle_style
     - yes
     - equilibrium angle for each type
   * - extract
     - all but test_fix_timestep
     - yes
     - list of keywords supported by the style's ``extract()`` method and their dimension
   * - init_stress
     - all but test_fix_timestep
     - yes
     - stress tensor after "run 0"
   * - init_forces
     - all but test_fix_timestep
     - yes
     - forces on atoms after "run 0"
   * - run_stress
     - all
     - no
     - stress tensor after the run (omitted by ``test_fix_timestep`` when the fix has no virial contribution)
   * - run_forces
     - all but test_fix_timestep
     - yes
     - forces on atoms after "run 4"
   * - run_pos
     - test_fix_timestep
     - yes
     - per-atom positions after the run
   * - run_vel
     - test_fix_timestep
     - yes
     - per-atom velocities after the run
   * - run_torque
     - test_fix_timestep
     - no
     - per-atom torques after the run (only when the atom style stores torque)
   * - global_scalar
     - test_fix_timestep
     - no
     - the global scalar output of the tested fix, if any
   * - global_vector
     - test_fix_timestep
     - no
     - the global vector output of the tested fix, if any

These reference files can be validated against the JSON schema file
``tools/json/force-style-test-schema.json`` with the ``check-jsonschema``
tool, which catches typos in keys, missing required keys, and values of the
wrong type.  For example, to validate all of them at once:

.. code-block:: sh

   check-jsonschema --schemafile tools/json/force-style-test-schema.json \
       unittest/force-styles/tests/*.yaml

See the :ref:`JSON support files <json>` section of the :doc:`Tools`
documentation for how to install ``check-jsonschema``.

The test program will read all this data from the YAML file and then
create a LAMMPS instance, apply the settings/commands from the YAML file
as needed and then issue a "run 0" command, write out a restart file, a
data file and a coeff file. The actual test will then compare computed
energies, stresses, and forces with the reference data, issue a "run 4"
command and compare to the second set of reference data.  This will be
run with both the newton_pair setting enabled and disabled and is
expected to generate the same results (allowing for some numerical
noise). Then it will restart from the previously generated restart and
compare with the reference and also start from the data file.  A final
check will use multi-cutoff r-RESPA (if supported by the pair style) at
a 1:1 split and compare to the Verlet results.  These sets of tests are
run with multiple test fixtures for accelerated styles: OPT, OPENMP and
INTEL (the latter two with 4 OpenMP threads enabled), and three mutually
exclusive KOKKOS fixtures selected by the active back end: the
``kokkos_omp`` fixture requires the KOKKOS package compiled with the
OpenMP back end and uses 4 OpenMP threads, while the ``kokkos_serial``
fixture only runs when the Serial back end is the sole back end of the
KOKKOS package (with any other back end enabled the host execution space
would not be Serial, so this configuration must be tested with a separate
build).  Both of these host fixtures skip when a GPU back end (CUDA, HIP,
SYCL) is enabled, since the KOKKOS package then must run on the GPU.  The
third fixture, ``kokkos_gpu``, is the complement: it runs only when a GPU
back end is enabled (using ``-k on g 1``) and is skipped on host-only
builds.  Because enabling the KOKKOS package with a GPU back end aborts
when no usable device is present, this fixture first probes for a
compatible GPU at runtime with ``Info::has_kokkos_gpu_device()`` (the
KOKKOS package analog of ``Info::has_gpu_device()`` for the GPU package)
and skips transparently when none is available, so the test suite can be
run unchanged on machines without a GPU.  For these tests the relative error
(epsilon) is lowered by a common factor due to the additional numerical
noise, but the tests are still comparing to the same reference data.

The KOKKOS fixtures also support the KOKKOS package compiled for reduced
precision with ``-D KOKKOS_PREC=mixed`` (compute in single precision,
accumulate in double precision) or ``-D KOKKOS_PREC=single``: the test
tolerance is then relaxed by a large additional factor, similar to what
is done for the mixed and single precision variants of the GPU package.
Individual tests can be skipped for a given fixture by listing the
fixture name in the ``skip_tests:`` field of the YAML file (e.g.
``skip_tests: kokkos_omp kokkos_serial kokkos_gpu``).  A skip entry may
also be qualified by the KOKKOS precision, e.g. ``kokkos_serial_single``
or ``kokkos_omp_mixed``, which skips the test only for that combination
of fixture and precision.  This is used for tests whose reference
quantities cannot be meaningfully compared in reduced precision, for
example global force totals that are the cancellation sum of large
per-atom contributions in a charge-neutral system.

The test fixture names accepted by ``skip_tests`` (each fixture runs the
corresponding variant or check and self-skips when its package or back end
is not available) are listed below.  Not every fixture exists for every
style kind (e.g. ``gpu``, ``intel``, and ``opt`` are used by the pair-style
tester, while ``numdiff`` is used by the bonded-style testers).

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Fixture
     - Description
   * - plain
     - the unmodified (base class) style, without a suffix
   * - omp
     - the ``/omp`` variant from the OPENMP package (run with 4 threads)
   * - intel
     - the ``/intel`` variant from the INTEL package
   * - opt
     - the ``/opt`` variant from the OPT package
   * - gpu
     - the ``/gpu`` variant from the GPU package
   * - kokkos_serial
     - the ``/kk`` variant from the KOKKOS package with a Serial-only build
   * - kokkos_omp
     - the ``/kk`` variant from the KOKKOS package with the OpenMP back end
       (run with 4 threads)
   * - kokkos_gpu
     - the ``/kk`` variant from the KOKKOS package with a GPU back end
   * - single
     - consistency check of the style's ``single()`` method against ``compute()``
   * - extract
     - check of the style's ``extract()`` keywords (base style)
   * - extract_omp
     - check of the style's ``extract()`` keywords (``/omp`` variant)
   * - numdiff
     - check of the forces against a numerical derivative of the energy

The ``kokkos_omp``, ``kokkos_serial``, and ``kokkos_gpu`` entries may be
qualified with ``_single`` or ``_mixed`` (e.g. ``kokkos_gpu_single``), as
noted above.

The ``tags:`` field of a YAML file lists keywords that classify a test or
request special handling from the test fixtures.  The fixtures query them
with ``TestConfig::has_tag()`` so that style-specific behavior is selected by
a descriptive tag instead of by hard-coded style names.  The recognized tags
are:

.. list-table::
   :header-rows: 1

   * - Tag
     - Purpose
   * - gpu_no_mixed
     - The GPU package variant of the style does not support mixed
       precision GPU mode; the ``gpu`` fixture skips the test when the
       GPU package was compiled for mixed precision.
   * - gpu_no_single
     - The GPU package variant of the style does not support single
       precision GPU mode (e.g. ``born/coul/long/cs/gpu``); the ``gpu``
       fixture skips the test when the GPU package was compiled for
       single precision.
   * - single_thread
     - The style cannot run correctly with more than one thread in the
       test (e.g.  ``dpd`` uses multiple per-thread pRNGs; ``snap`` and
       ``pace`` due to their implementation), so the threaded fixtures
       (``omp``, ``intel``, ``kokkos_omp``) run it with a single thread.
   * - no_respa
     - The ``fix_timestep`` tester does not exercise this style under
       :doc:`run_style respa <run_style>`: rigid fixes need additional
       work to test correctly with r-RESPA, ``fix nve/limit`` and ``fix
       recenter`` do not support it, stochastic integrators and barostats
       (``brownian``, ``gjf``, ``press/langevin``) draw their random
       numbers differently under r-RESPA, and velocity-dependent forcing
       fixes (``viscous``, ``accelerate/cos``) and the isokinetic ``nvk``
       integrator follow a different trajectory under r-RESPA - in all of
       these cases the verlet and r-RESPA runs cannot match.
   * - no_reset_dt
     - The ``fix_timestep`` tester does not exercise a timestep change
       for this style.  The fix rejects a timestep reset (its
       ``Fix::reset_dt()`` raises an error, e.g. :doc:`fix move
       <fix_move>`), which would otherwise abort the test.
   * - ellipsoid
     - The test includes ellipsoids and thus requires :doc:`fix
       nve/asphere <fix_nve_asphere>`.
   * - spica_pair
     - The test setup uses ``pair_style lj/spica`` instead of the
       default ``pair_style zero`` (required by the ``spica`` angle
       style).
   * - slow
     - The test runs significantly longer than others and ``ctest -LE
       slow`` would skip it.
   * - noWindows
     - Indicates that this test must be skipped on Windows; use
       ``ctest -LE noWindows``
   * - unstable
     - The test exhibits numerically unstable behavior on some
       platforms, e.g. ARM64; Until a proper correction is found, tests
       can be skipped with ``ctest -LE unstable``.
   * - generated
     - Indicates that a test input was regenerated. *Remove* after
       confirming the correctness of the updated YAML file.

Additional tests will check whether all listed extract keywords are
supported and have the correct dimensionality and the final set of tests
will set up a few pairs of atoms explicitly and in such a fashion that
the forces on the atoms computed from ``Pair::compute()`` will match
individually with the results from ``Pair::single()``, if the pair style
does support that functionality.

With this scheme a large fraction of the code of any tested pair style
will be executed and consistent results are required for different
settings and between different accelerated pair style variants and the
base class, as well as for computing individual pairs through the
``Pair::single()`` method where supported.

The ``test_pair_style`` tester is used with 4 categories of test inputs:

- pair styles compatible with molecular systems using bonded
  interactions and exclusions.  For pair styles requiring a KSpace style
  the KSpace computations are disabled.  The YAML files match the
  pattern "mol-pair-\*.yaml" and the tests are correspondingly labeled
  with "MolPairStyle:\*"
- pair styles not compatible with the previous input template.
  The YAML files match the pattern "atomic-pair-\*.yaml" and the tests are
  correspondingly labeled with "AtomicPairStyle:\*"
- manybody pair styles.
  The YAML files match the pattern "atomic-pair-\*.yaml" and the tests are
  correspondingly labeled with "AtomicPairStyle:\*"
- kspace styles.
  The YAML files match the pattern "kspace-\*.yaml" and the tests are
  correspondingly labeled with "KSpaceStyle:\*". In these cases a compatible
  pair style is defined, but the computation of the pair style contributions
  is disabled.

The ``test_bond_style``, ``test_angle_style``, ``test_dihedral_style``, and
``test_improper_style`` tester programs are set up in a similar fashion and
share support functions with the pair style tester.  The final group of
tests in this section is for fix styles that add/manipulate forces and
velocities, e.g. for time integration, thermostats and more.

Adding a new test is easiest done by copying and modifying an existing YAML
file for a style that is similar to one to be tested.  The file name should
follow the naming conventions described above and after copying the file,
the first step is to replace the style names where needed.  The coefficient
values do not have to be meaningful, just in a reasonable range for the
given system.  It does not matter if some forces are large, for as long as
they do not diverge.

The template input files define a large number of index variables at the top
that can be modified inside the YAML file to control the behavior.  For example,
if a pair style requires a "newton on" setting, the following can be used in
as the "pre_commands" section:

.. code-block:: yaml

   pre_commands: ! |
     variable newton_pair delete
     variable newton_pair index on

And for a pair style requiring a kspace solver the following would be used as
the "post_commands" section:

.. code-block:: yaml

   post_commands: ! |
     pair_modify table 0
     kspace_style pppm/tip4p 1.0e-6
     kspace_modify gewald 0.3
     kspace_modify compute no

Note that this disables computing the kspace contribution, but still will run
the setup.  The "gewald" parameter should be set explicitly to speed up the run.
For styles with long-range electrostatics, typically two tests are added one using
the (slower) analytic approximation of the erfc() function and the other using
the tabulated coulomb, to test both code paths.  The reference results in the YAML
files then should be compared manually, if they agree well enough within the limits
of those two approximations.

The ``test_pair_style`` and equivalent programs have special command-line options
to update the YAML files. Running a command like

.. code-block:: bash

   test_pair_style mol-pair-lennard_mdf.yaml -g new.yaml

will read the settings from the ``mol-pair-lennard_mdf.yaml`` file and then compute
the reference data and write a new file with to ``new.yaml``.  If this step fails,
there are likely some (LAMMPS or YAML) syntax issues in the YAML file that need to
be resolved and then one can compare the two files to see if the output is as expected.

It is also possible to do an update in place with:

.. code-block:: bash

   test_pair_style mol-pair-lennard_mdf.yaml -u

And one can finally run the full set of tests with:

.. code-block:: bash

   test_pair_style mol-pair-lennard_mdf.yaml

This will just print a summary of the groups of tests.  When using the "-v" flag
the test will also keep any LAMMPS output and when using the "-s" flag, there
will be some statistics reported on the relative errors for the individual checks
which can help to figure out what would be a good choice of the epsilon parameter.
It should be as small as possible to catch any unintended side effects from changes
elsewhere, but large enough to accommodate the numerical noise due to the implementation
of the potentials and differences in compilers.

.. note::

   These kinds of tests can be very sensitive to compiler optimization and
   thus the expectation is that they pass with compiler optimization turned
   off. When compiler optimization is enabled, there may be some failures, but
   one has to carefully check whether those are acceptable due to the enhanced
   numerical noise from reordering floating-point math operations or due to
   the compiler mis-compiling the code. That is not always obvious.


Tests for granular (DEM) models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionadded:: 4Jul2026

The ``unittest/granular`` folder contains a YAML-driven test suite for
discrete element method (DEM) / granular models, built in the same spirit
as the force-style tests above but specialized for time-resolved
trajectories of small granular systems.

Currently, there are 11 test programs. This set of unit tests is still a
work-in-progress and the tests have not yet been thoroughly vetted. Tests
may be added, updated, or removed. The first six test programs,
``test_dem_01`` through ``test_dem_06``, reproduce the test surface of the
MFiX-DEM verification studies of :ref:`Garg et al. <dem_Garg2012>` (the
individual cases are also described in the `MFiX-DEM manual
<https://mfix.netl.doe.gov/doc/vvuq-manual/main/html/dem/index.html>`_).
``test_dem_07`` through ``test_dem_11`` add coverage of benchmark cases
from the granular literature: rolling resistance, cohesion, and
two-particle collisions follow the software-agnostic DEM benchmark of
:ref:`Mohajeri et al. <dem_Mohajeri2024>`, and the bulk angle of repose
and multi-sphere clump cases follow the round-robin study of
:ref:`Saomoto et al. <dem_Saomoto2023>`.  The particle-impact-level cases
-- oblique wall impact (``test_dem_05``), two-sphere and spinning-sphere
collisions (``test_dem_09``) and the elastic Hertzian normal impact
(``test_dem_11``) -- follow the benchmark of :ref:`Chung and Ooi
<dem_Chung2011>`.  The test programs are:

.. list-table::
   :header-rows: 1

   * - Program
     - Scenario
   * - ``test_dem_01``
     - a freely falling particle bouncing off a wall
   * - ``test_dem_02``
     - a particle bouncing repeatedly (convergence to the hard-sphere limit)
   * - ``test_dem_03``
     - two stacked particles in continuous compression between two walls
   * - ``test_dem_04``
     - a sphere sliding then rolling without slipping on a rough surface
   * - ``test_dem_05``
     - an oblique collision of a sphere (and a superellipsoid) with a wall
   * - ``test_dem_06``
     - a single particle settling to its terminal velocity under fluid drag
   * - ``test_dem_07``
     - a spinning sphere damped to rest by rolling resistance (``rolling sds``)
   * - ``test_dem_08``
     - cohesive/adhesive contact: the DMT and JKR pull-off force
   * - ``test_dem_09``
     - two-sphere head-on, oblique (shear), and spinning-sphere collisions
   * - ``test_dem_10``
     - bulk behavior: a settling pile, the angle of repose, and a rigid clump
   * - ``test_dem_11``
     - elastic Hertzian normal impact (peak contact mechanics)

Every test program shares the same driver logic, implemented in
``unittest/granular/test_dem_common.cpp`` and compiled into the
``granular_tests`` support library; each ``test_dem_0N.cpp`` only contains
the two GoogleTest fixtures (``newton_on`` and ``newton_off``).  As with the
force-style tests, the reference systems are defined by YAML files in the
``unittest/granular/tests`` folder and registered as CTest cases by their
file name (``dem0N-*.yaml`` becomes test ``DEM0N:*``); adding or removing a
YAML file requires re-running CMake.

Unlike the force-style tests, the entire system is built *from the YAML
file* rather than from a fixed input template.  A YAML file provides an
optional ``variables`` block (emitted as :doc:`index variables <variable>`
so they can be substituted as ``${name}`` anywhere in the command strings),
``pre_commands`` that create the geometry, ``pair_style`` / ``pair_coeff``
that select the contact model, and ``post_commands`` that add the
integrator, gravity, walls and drag.  The trajectory is then advanced in a
sequence of ``run_segments`` and, after each segment, the per-atom
positions, velocities, torques, angular velocities (spheres) and angular
momenta (ellipsoids/superellipsoids) are compared against the recorded
reference.  A minimal example (``dem01-hooke-3d-si.yaml``) looks like:

.. code-block:: yaml

   ---
   lammps_version: 30 Mar 2026
   tags: granular
   epsilon: 1e-10
   prerequisites: ! |
     atom sphere
     pair gran/hooke
   variables: ! |
     knorm 1.0e4
     gnorm 10.0
     diam 0.2
     dens 2600.0
     grav 9.81
     z0 0.5
   pre_commands: ! |
     units si
     dimension 3
     boundary f f f
     atom_style sphere
     region box block -0.5 0.5 -0.5 0.5 0.0 1.0 units box
     create_box 1 box
     create_atoms 1 single 0.0 0.0 ${z0} units box
     set group all diameter ${diam} density ${dens}
     comm_modify vel yes
     timestep 0.001
   pair_style: gran/hooke ${knorm} NULL ${gnorm} NULL 0.0 0
   pair_coeff: ! |
     * *
   post_commands: ! |
     fix grav all gravity ${grav} vector 0.0 0.0 -1.0
     fix integr all nve/sphere
     fix zwall all wall/gran hooke ${knorm} NULL ${gnorm} NULL 0.0 0 zplane 0.0 NULL
   run_segments: ! |
     250 150 300
   analytic_enable: yes
   analytic_model: freefall
   analytic_tol: 1.0e-9
   analytic_segment: 0
   # run_pos / run_vel / run_torque / run_omega / run_angmom blocks follow

The following table describes the available keys:

.. list-table::
   :header-rows: 1

   * - Key:
     - Description:
   * - epsilon
     - relative precision required for the recorded (regression) reference data
   * - prerequisites
     - list of style kind / style name pairs required to run the test
   * - variables
     - name/value pairs exposed as ``${name}`` index variables for substitution
   * - pre_commands
     - commands that build the geometry (units, box, atoms, ``set``, timestep)
   * - pair_style / pair_coeff
     - the particle-particle contact model
   * - post_commands
     - fixes added after the geometry (integrator, gravity, walls, drag)
   * - run_segments
     - whitespace-separated list of run lengths; state is captured after each
   * - run_pos, run_vel
     - reference positions and velocities, as ``segment tag x y z`` rows
   * - run_torque, run_omega, run_angmom
     - reference torque / angular velocity / angular momentum (when applicable)
   * - analytic_enable
     - ``yes`` to also assert a closed-form (analytic) model
   * - analytic_model
     - which analytic model to evaluate (see below)
   * - analytic_tol
     - relative tolerance for the analytic assertion (looser than ``epsilon``)
   * - analytic_segment
     - run segment at which the analytic model is checked (``-1`` means the last)
   * - analytic_only
     - ``yes`` to record/check *only* the analytic model and skip the per-atom
       regression (for chaotic bulk tests; see below)

The per-atom reference blocks use a ``segment tag x y z`` row format, so a
single block holds the data for all run segments and the row order does not
matter.  Because granular/atomic systems do not build an atom map by
default, the reference generator iterates over local atoms by tag rather
than calling ``Atom::map()``.

Each test runs as a pure regression check (the recorded data is reproduced
to within ``epsilon``) under both the ``newton on`` and ``newton off``
fixtures, which are expected to give identical results.  In addition, a
test may opt in to an *analytic* check that compares a derived quantity
against a closed-form solution implemented in
``unittest/granular/test_analytic_models.cpp``.  The analytic tolerance is
deliberately loose, because the soft-sphere DEM result only approaches the
idealized (hard-sphere or instantaneous-contact) solution.  The models
currently implemented are:

.. list-table::
   :header-rows: 1

   * - Model:
     - Checks:
   * - freefall
     - ballistic motion before contact: :math:`z = z_0 - g t^2/2`, :math:`v_z = -g t`
   * - bounce_height
     - hard-sphere apex after the k-th bounce :math:`h_k = r + e^{2k}(h_0 - r)`
   * - stack_energy
     - conservation of total mechanical energy for an elastic two-particle stack
   * - slip_cessation
     - rolling-without-slipping limit :math:`u = 5 u_0/7`, :math:`\omega = u/r`
   * - oblique_impact
     - gross-sliding rebound :math:`v_x' = v_x - \mu(1+e)v_z`, :math:`\omega_y = \tfrac{5}{2}\mu(1+e)v_z/r`
   * - terminal_velocity_linear
     - Stokes drag terminal velocity :math:`v_{term} = m g/\gamma`
   * - terminal_velocity_schiller_naumann
     - Schiller-Naumann terminal velocity from :math:`m g = \tfrac{1}{2} C_d \rho_g \pi r^2 v^2`
   * - rolling_decay
     - linear spin-down under rolling resistance: :math:`\omega = \omega_0 - \tfrac{5 \mu_r g}{2 r} t`
   * - pulloff_dmt
     - DMT pull-off force at contact :math:`|F| = 4 \pi \gamma R_{\mathrm{eff}}`
   * - collision_restitution
     - two-sphere momentum conservation and restitution :math:`e = -(v_1'-v_2')/(v_1-v_2)`
   * - angle\_of\_repose
     - measured heap slope :math:`\arctan(z_{\max}/r_{\max})` lies within a ``[lo, hi]`` band
   * - hertz\_normal\_impact
     - Hertzian peak energy balance :math:`\tfrac{1}{2}\mu_{red} V_{rela}^2 = \tfrac{2}{5} P_{max}\alpha_{max}`
   * - spin\_impact
     - gross-sliding rebound of a spinning sphere: :math:`v_x' = \mu(1+e)v_n`, :math:`\omega_y' = \omega_0 - \tfrac{5}{2}\mu(1+e)v_n/r`
   * - spin\_no\_friction
     - counter-spinning spheres with zero contact slip keep their spin and gain no tangential velocity

``test_dem_06`` exercises both :doc:`fix viscous <fix_viscous>` (linear
Stokes drag) and the :doc:`fix viscous/nonlinear <fix_viscous_nonlinear>`
style that was added together with these tests for the Schiller-Naumann
drag correlation.

Analytic-only (chaotic bulk)
""""""""""""""""""""""""""""

Most tests are bit-for-bit
regressions that reproduce identically under ``newton on`` and ``newton
off``.  A few bulk scenarios -- notably the angle-of-repose pile in
``test_dem_10`` -- are *chaotic*: a long pour-and-settle trajectory amplifies
the round-off differences between summation orders, so the per-atom state is
not reproducible across ``newton`` settings or platforms even though the bulk
observable (the heap angle) is robust.  Such a YAML sets ``analytic_only:
yes``. The generator then records no per-atom reference blocks, and the
driver checks only the analytic model.  ``test_dem_10`` pairs this with a
short, deterministic ``dem10-settle-*`` regression (a small lattice block
relaxing into contact) and a deterministic ``dem10-clump-*`` case (a
:doc:`fix rigid/small <fix_rigid>` tetrahedral clump bouncing on a granular
wall) so the bit-for-bit code path is still covered.

Adding a new reference (YAML) file
""""""""""""""""""""""""""""""""""

Copy an existing ``dem0N-*.yaml``
for a similar scenario, adjust the ``variables``, ``pre_commands``,
``pair_style``/``pair_coeff`` and ``post_commands`` for the new model, and
give it a new name matching the ``dem0N-*.yaml`` pattern of the test program
it belongs to.  Leave out the reference data blocks initially, then
(re)generate them in place with:

.. code-block:: bash

   TEST_ARGS=-u ctest -R DEM0N:myvariant

or by running the driver directly (``test_dem_0N dem0N-myvariant.yaml -u``).
Do **not** write the generated file to a sibling ``dem0N-*.yaml`` name (for
example with the ``-g newfile.yaml`` option pointing into the ``tests``
folder), because the ``CONFIGURE_DEPENDS`` glob would then register it as an
extra, stale test.  After adding the file, re-run CMake so the new test is
registered, then verify it with ``ctest -V -R DEM0N:myvariant`` (the ``-s``
option of the driver reports per-quantity error statistics, which helps when
choosing ``epsilon`` and the analytic tolerance).

Adding a new test program
"""""""""""""""""""""""""

Create ``test_dem_0N.cpp`` as a thin copy of an existing one (only the
GoogleTest suite name changes), add an
``add_executable``/``register_dem_tests`` pair to
``unittest/granular/CMakeLists.txt``, and add the corresponding
``dem0N-*.yaml`` reference files.  If the new scenario needs a
closed-form check, add a named model to ``test_analytic_models.cpp``
that reads its parameters from the ``variables`` block (and reads
masses, radii, etc. from the live LAMMPS instance to avoid depending on
derived quantities) and assert it with ``EXPECT_LE`` on the relative
error.

References
""""""""""

.. _dem_Garg2012:

**(Garg et al., 2012)** R. Garg, J. Galvin, T. Li, and S. Pannala,
Open-source MFIX-DEM software for gas-solids flows: Part I -- Verification
studies, Powder Technology, 220, 122-137 (2012),
https://doi.org/10.1016/j.powtec.2011.09.019

.. _dem_Mohajeri2024:

**(Mohajeri et al., 2024)** M. J. Mohajeri, C. Coetzee, and D. L. Schott,
A software-agnostic benchmark for DEM simulation of cohesive and
non-cohesive materials, Powder Technology, 447, 120136 (2024),
https://doi.org/10.1016/j.powtec.2024.120136

.. _dem_Saomoto2023:

**(Saomoto et al., 2023)** H. Saomoto, N. Kikkawa, S. Moriguchi, Y. Nakata,
et al., Round robin test on angle of repose: DEM simulation results
collected from 16 groups around the world, Soils and Foundations, 63,
101272 (2023), https://doi.org/10.1016/j.sandf.2023.101272

.. _dem_Chung2011:

**(Chung and Ooi, 2011)** Y. C. Chung and J. Y. Ooi, Benchmark tests for
verifying discrete element modelling codes at particle impact level,
Granular Matter, 13, 643-656 (2011),
https://doi.org/10.1007/s10035-011-0277-0


Tests for programs in the tools folder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``unittest/tools`` folder contains tests for programs in the
``tools`` folder.  This currently only contains tests for the LAMMPS
shell, which are implemented as a python scripts using the ``unittest``
Python module and launching the tool commands through the ``subprocess``
Python module.


Troubleshooting failed unit tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are by default no unit tests for newly added features (e.g. pair,
fix, or compute styles) unless your pull request also includes tests for
these added features.  If you are modifying some existing LAMMPS
features, you may see failures for existing tests, if your modifications
have some unexpected side effects or your changes render the existing
test invalid.  If you are adding an accelerated version of an existing
style, then only tests for INTEL, KOKKOS (with the OpenMP or Serial host
back ends, depending on how the KOKKOS package was configured), OPENMP,
and OPT will be run automatically.  Tests for the GPU package are time
consuming and thus are only run *after* a merge, or when a special
label, ``gpu_unit_tests`` is added to the pull request.  After the test
has started, it is often best to remove the label since every PR
activity will re-trigger the test (that is a limitation of triggering a
test with a label).  Support for unit tests using KOKKOS with GPU
acceleration is currently not supported.

When you see a failed build on GitHub, click on ``Details`` to be taken
to the corresponding LAMMPS Jenkins CI web page.  Click on the "Exit"
symbol near the ``Logout`` button on the top right of that page to go to
the "classic view".  In the classic view, there is a list of the
individual runs that make up this test run (they are shown but cannot be
inspected in the default view).  You can click on any of those.
Clicking on ``Test Result`` will display the list of failed tests. Click
on the "Status" column to sort the tests based on their Failed or Passed
status.  Then click on the failed test to expand its output.

For example, the following output snippet shows the failed unit test

.. code-block:: console

   [ RUN      ] PairStyle.gpu
   /home/builder/workspace/dev/pull_requests/ubuntu_gpu/unit_tests/cmake_gpu_opencl_mixed_smallbig_clang_static/unittest/force-styles/test_main.cpp:63: Failure
   Expected: (err) <= (epsilon)
   Actual: 0.00018957912910606503 vs 0.0001
   Google Test trace:
   /home/builder/workspace/dev/pull_requests/ubuntu_gpu/unit_tests/cmake_gpu_opencl_mixed_smallbig_clang_static/unittest/force-styles/test_main.cpp:56: EXPECT_FORCES: init_forces (newton off)
   /home/builder/workspace/dev/pull_requests/ubuntu_gpu/unit_tests/cmake_gpu_opencl_mixed_smallbig_clang_static/unittest/force-styles/test_main.cpp:64: Failure
   Expected: (err) <= (epsilon)
   Actual: 0.00022892713393549854 vs 0.0001

The failed assertions provide line numbers in the test source
(e.g. ``test_main.cpp:56``), from which one can understand what
specific assertion failed.

Note that the force style engine runs one of a small number of systems
in a rather off-equilibrium configuration with a few atoms for a few
steps, writes data and restart files, uses :doc:`the clear command
<clear>` to reset LAMMPS, and then runs from those files with different
settings (e.g. newton on/off) and integrators (e.g. verlet vs. respa).
Beyond potential issues/bugs in the source code, the mismatch between
the expected and actual values could be that force arrays are not
properly cleared between multiple run commands or that class members are
not correctly initialized or written to or read from a data or restart
file.

While the epsilon (relative precision) for a single, `IEEE 754 compliant
<https://en.wikipedia.org/wiki/IEEE_754>`_, double precision floating
point operation is at about 2.2e-16, the achievable precision for the
tests is lower due to most numbers being sums over intermediate results
for which the non-associativity of floating point math leads to larger
errors.  As a rule of thumb, the test epsilon can often be in the range
5.0e-14 to 1.0e-13.  But for "noisy" force kernels, e.g. those a larger
amount of arithmetic operations involving `exp()`, `log()` or `sin()`
functions, and also due to the effect of compiler optimization or differences
between compilers or platforms, epsilon may need to be further relaxed,
sometimes epsilon can be relaxed to 1.0e-12. If interpolation or lookup
tables are used, epsilon may need to be set to 1.0e-10 or even higher.
For tests of accelerated styles, the per-test epsilon is multiplied
by empirical factors that take into account the differences in the order
of floating point operations or that some or most intermediate operations
may be done using approximations or with single precision floating point
math.

To rerun a failed unit test individually, change to the ``build`` directory
and run the test with verbose output. For example,

.. code-block:: bash

    env TEST_ARGS=-v ctest -R ^MolPairStyle:lj_cut_coul_long -V

``ctest`` with the ``-V`` flag also shows the exact command
of the test. One can then use ``gdb --args`` to further debug and
catch exceptions with the test command, for example,

.. code-block:: bash

    gdb --args /path/to/lammps/build/test_pair_style /path/to/lammps/unittest/force-styles/tests/mol-pair-lj_cut_coul_long.yaml


It is recommended to configure the build with ``-D
BUILD_SHARED_LIBS=on`` and use a custom linker to shorten the build time
during recompilation.  Installing `ccache` in your development
environment helps speed up recompilation by caching previous
compilations and detecting when the same compilation is being done
again.  Please see :doc:`Build_development` for further details.
