Adding tests for unit testing
-----------------------------

This section discusses adding or expanding tests for the unit test
infrastructure included into the LAMMPS source code distribution.
Unlike example inputs, unit tests focus on testing the "local" behavior
of individual features, tend to run very fast, and should be set up to
cover as much of the added code as possible.  When contributing code to
the distribution, the LAMMPS developers will appreciate if additions
to the integrated unit test facility are included.

Given the complex nature of MD simulations where many operations can
only be performed when suitable "real" simulation environment has been
set up, not all tests will be unit tests in the strict definition of
the term.  They are rather be executed on a more abstract level by issuing
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
``test_tokenizer.cpp`` and links it with the googletest libraries and the
LAMMPS library as well as it uses the ``main()`` function from the
GMock library of googletest.  The third line registers the executable
as a test program to be run from ``ctest`` under the name ``Tokenizer``.

The test executable itself will execute multiple individual tests
through the googletest framework. In this case each test consists of
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
and does not use the googletest framework and thus not representative.
The other test sources, however, can serve as guiding examples for
additional tests.

Tests for individual LAMMPS commands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These tests are a bit more complex as they require to first create a
:cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class instance and then use the
:doc:`C++ API <Cplusplus>` to pass individual commands to that LAMMPS
instance.  For that reason these tests use a googletest "test fixture",
i.e. a class derived from ``testing::Test`` that will create (and
delete) the required LAMMPS class instance for each set of tests in a
TEST_F() function.  Please see the individual source files for different
examples of setting up suitable test fixtures.  Here is an example for
implementing a test using a fixture by first checking the default
value and then issuing LAMMPS commands and checking whether they
have the desired effect:

.. code-block:: c++

   TEST_F(SimpleCommandsTest, ResetTimestep)
   {
       ASSERT_EQ(lmp->update->ntimestep, 0);

       if (!verbose) ::testing::internal::CaptureStdout();
       lmp->input->one("reset_timestep 10");
       if (!verbose) ::testing::internal::GetCapturedStdout();
       ASSERT_EQ(lmp->update->ntimestep, 10);

       if (!verbose) ::testing::internal::CaptureStdout();
       lmp->input->one("reset_timestep 0");
       if (!verbose) ::testing::internal::GetCapturedStdout();
       ASSERT_EQ(lmp->update->ntimestep, 0);
   }

Please note the use of the (global) verbose variable to control whether
the LAMMPS command will be silent by capturing the output or not.  In
the default case, verbose == false, the test output will be compact and
not mixed with LAMMPS output. However setting the verbose flag (via
setting the TEST_ARGS environment variable, ``TEST_ARGS=-v``) can be
helpful to understand why tests fail unexpectedly.
   
Another complexity of these tests stems from the need to capture
situations where LAMMPS will stop with an error, i.e. handle so-called
"death tests".  Here the LAMMPS code will operate differently depending
on whether it was configured to throw C++ exceptions on errors or call
either ``exit()`` or ``MPI_Abort()``.  In the latter case, the test code
also needs to detect whether LAMMPS was compiled with the OpenMPI
library, as OpenMPI is **only** compatible the death test options of the
googletest library when C++ exceptions are enabled; otherwise those
"death tests" must be skipped to avoid reporting bogus failures.  The
specifics of this step are implemented in the TEST_FAILURE()
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
   * - ``test_kim_commands.cpp``
     - KimCommands
     - Tests for several commands from the :ref:`KIM package <PKG-KIM>`
   * - ``test_reset_ids.cpp``
     - ResetIDs
     - Tests to validate the :doc:`reset_atom_ids <reset_atom_ids>` and :doc:`reset_mol_ids <reset_mol_ids>` commands


Adding tests for the C-style library interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for the Python module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for the Fortran interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for the C++-style library interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for styles computing or modifying forces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


Adding tests for programs in the tools folder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
