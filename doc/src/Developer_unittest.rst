Adding tests for unit testing
-----------------------------

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
   :widths: auto
   :align: left

   * - File name:
     - Test name:
     - Description:
   * - ``test_argutils.cpp``
     - ArgInfo
     - Tests for ``ArgInfo`` class used by LAMMPS
   * - ``test_fmtlib.cpp``
     - FmtLib
     - Tests for ``fmtlib::`` functions used by LAMMPS
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

.. code-block::

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

Another complexity of these tests stems from the need to capture
situations where LAMMPS will stop with an error, i.e. handle so-called
"death tests".  Here the LAMMPS code will operate differently depending
on whether it was configured to throw C++ exceptions on errors or call
either ``exit()`` or ``MPI_Abort()``.  In the latter case, the test code
also needs to detect whether LAMMPS was compiled with the OpenMPI
library, as OpenMPI is **only** compatible the death test options of the
GoogleTest library when C++ exceptions are enabled; otherwise those
"death tests" must be skipped to avoid reporting bogus failures.  The
specifics of this step are implemented in the ``TEST_FAILURE()``
macro. These tests operate by capturing the screen output when executing
the failing command and then comparing that with a provided regular
expression string pattern.  Example:

.. code-block:: C++

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
   * - ``test_reset_ids.cpp``
     - ResetIDs
     - Tests to validate the :doc:`reset_atom_ids <reset_atom_ids>` and :doc:`reset_mol_ids <reset_mol_ids>` commands


Tests for the C-style library interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tests for validating the LAMMPS C-style library interface are in the
``unittest/c-library`` folder.  They are implemented in either way used
for utility functions and for LAMMPS commands, but use the functions
implemented in the ``src/library.cpp`` file as much as possible.  There
may be some overlap with other tests, but only in as much as is required
to test the C-style library API.  The tests are distributed over
multiple test programs which tries to match the grouping of the
functions in the source code and :ref:`in the manual <lammps_c_api>`.

This group of tests also includes tests invoking LAMMPS in parallel
through the library interface, provided that LAMMPS was compiled with
MPI support.  These include tests where LAMMPS is run in multi-partition
mode or only on a subset of the MPI world communicator.  The CMake
script code for adding this kind of test looks like this:

.. code-block:: CMake

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
to also implement test wrappers in C++ that will call fortran functions
which provide a C function interface through ISO_C_BINDINGS that will in
turn call the functions in the LAMMPS Fortran module.  This part of the
unit tests is incomplete since the Fortran module it is based on is
incomplete as well.

Tests for the C++-style library interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tests in the ``unittest/cplusplus`` folder are somewhat similar to
the tests for the C-style library interface, but do not need to test the
several convenience and utility functions that are only available through
the C-style interface.  Instead it can focus on the more generic features
that are used internally.  This part of the unit tests is currently still
mostly in the planning stage.

Tests for reading and writing file formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``unittest/formats`` folder contains test programs for reading and
writing files like data files, restart files, potential files or dump files.
This covers simple things like the file i/o convenience functions in the
``utils::`` namespace to complex tests of atom styles where creating and
deleting of atoms with different properties is tested in different ways
and through script commands or reading and writing of data or restart files.

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

The following table describes the available keys and their purpose for
testing pair styles:

.. list-table::
   :header-rows: 1

   * - Key:
     - Description:
   * - lammps_version
     - LAMMPS version used to last update the reference data
   * - date_generated
     - date when the file was last updated
   * - epsilon
     - base value for the relative precision required for tests to pass
   * - prerequisites
     - list of style kind / style name pairs required to run the test
   * - pre_commands
     - LAMMPS commands to be executed before the input template file is read
   * - post_commands
     - LAMMPS commands to be executed right before the actual tests
   * - input_file
     - LAMMPS input file template based on pair style zero
   * - pair_style
     - arguments to the pair_style command to be tested
   * - pair_coeff
     - list of pair_coeff arguments to set parameters for the input template
   * - extract
     - list of keywords supported by ``Pair::extract()`` and their dimension
   * - natoms
     - number of atoms in the input file template
   * - init_vdwl
     - non-Coulomb pair energy after "run 0"
   * - init_coul
     - Coulomb pair energy after "run 0"
   * - init_stress
     - stress tensor after "run 0"
   * - init_forces
     - forces on atoms after "run 0"
   * - run_vdwl
     - non-Coulomb pair energy after "run 4"
   * - run_coul
     - Coulomb pair energy after "run 4"
   * - run_stress
     - stress tensor after "run 4"
   * - run_forces
     - forces on atoms after "run 4"

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
run with multiple test fixtures for accelerated styles (OPT, OPENMP,
INTEL) and for the latter two with 4 OpenMP threads enabled.  For
these tests the relative error (epsilon) is lowered by a common factor
due to the additional numerical noise, but the tests are still comparing
to the same reference data.

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
``Pair::single()`` where supported.

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

The ``test_bond_style`` and ``test_angle_style`` are set up in a similar
fashion and share support functions with the pair style tester.  The final
group of tests in this section is for fix styles that add/manipulate forces
and velocities, e.g. for time integration, thermostats and more.

Adding a new test is easiest done by copying and modifying an existing test
for a style that is similar to one to be tested.  The file name should follow
the naming conventions described above and after copying the file, the first
step is to replace the style names where needed.  The coefficient values
do not have to be meaningful, just in a reasonable range for the given system.
It does not matter if some forces are large, for as long as they do not diverge.

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
the tabulated coulomb, to test both code paths. The reference results in the YAML
files then should be compared manually, if they agree well enough within the limits
of those two approximations.

The ``test_pair_style`` and equivalent programs have special command line options
to update the YAML files. Running a command like

.. code-block:: bash

   $ test_pair_style mol-pair-lennard_mdf.yaml -g new.yaml

will read the settings from the ``mol-pair-lennard_mdf.yaml`` file and then compute
the reference data and write a new file with to ``new.yaml``.  If this step fails,
there are likely some (LAMMPS or YAML) syntax issues in the YAML file that need to
be resolved and then one can compare the two files to see if the output is as expected.

It is also possible to do an update in place with:

.. code-block:: bash

   $ test_pair_style mol-pair-lennard_mdf.yaml -u

And one can finally run the full set of tests with:

.. code-block:: bash

   $ test_pair_style mol-pair-lennard_mdf.yaml

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


Tests for programs in the tools folder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``unittest/tools`` folder contains tests for programs in the
``tools`` folder.  This currently only contains tests for the LAMMPS
shell, which are implemented as a python scripts using the ``unittest``
Python module and launching the tool commands through the ``subprocess``
Python module.
