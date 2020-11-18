Adding tests for unit testing
-----------------------------

This section discusses adding or expanding tests for the unit test
infrastructure included into the LAMMPS source code distribution.
Unlike example inputs, unit tests focus on testing the "local" behavior
of individual features, tend to run very fast, and should be set up to
cover as much of the added code as possible.  When contributing code to
the distribution, the LAMMPS developers will appreciate if additions
to the integrated unit test facility are included.

Given the complex nature of MD simulations (i.e. many operations can
only be performed when suitable "real" simulation environment has been
set up) not all tests will be unit tests in the strict definition of
the term, but rather be executed on a more abstract level by issuing
LAMMPS script commands and then inspecting the changes to the internal
data.  For some classes of tests, generic test programs have been
written that can be applied to parts of LAMMPS that use the same
interface (via polymorphism) and those are driven by input files, so
tests can be added by simply adding more of those input files.


Adding tests for styles computing forces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for individual LAMMPS commands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for the C-style library interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for the Python module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for the Fortran interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for the C++-style library interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for utility functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding tests for programs in the tools folder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
