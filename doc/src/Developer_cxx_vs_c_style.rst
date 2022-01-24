Code design
-----------

This section discusses some of the code design choices in LAMMPS and
overall strategy in order to assist developers to write new code that
will fit well with the remaining code.  Please see the section on
:doc:`Requirements for contributed code <Modify_style>` for more
specific recommendations and guidelines.  Here the focus is on overall
strategy and discussion of some relevant C++ programming language
constructs.

Historically, the basic design philosophy of the LAMMPS C++ code was
that of a "C with classes" style.  The was motivated by the desire to
make it easier to modify LAMMPS for people without significant training
in C++ programming and by trying to use data structures and code constructs
that somewhat resemble the previous implementation(s) in Fortran.
A contributing factor for this choice also was that at the time the
implementation of C++ compilers was not always very mature and some of
the advanced features contained bugs or were not functioning exactly
as the standard required; plus there was some disagreement between
compiler vendors about how to interpret the C++ standard documents.

However, C++ compilers have advanced a lot since then and with the
transition to requiring the C++11 standard in 2020 as the minimum C++ language
standard for LAMMPS, the decision was made to also replace some of the
C-style constructs with equivalent C++ functionality, either from the
C++ standard library or as custom classes or function, in order to
improve readability of the code and to increase code reuse through
abstraction of commonly used functionality.


Object oriented code
^^^^^^^^^^^^^^^^^^^^

LAMMPS is designed to be an object oriented code, that is each simulation
is represented by an instance of the LAMMPS class.  When running in parallel,
of course, each MPI process will create such an instance.  This can be seen
in the ``main.cpp`` file where the core steps of running a LAMMPS simulation
are the following 3 lines of code:

.. code-block:: C++

    LAMMPS *lammps = new LAMMPS(argc, argv, lammps_comm);
    lammps->input->file();
    delete lammps;

The first line creates a LAMMPS class instance and passes the command line arguments
and the global communicator to its constructor.  The second line tells the LAMMPS
instance to process the input (either from standard input or the provided input file)
until the end.  And the third line deletes that instance again.  The remainder of
the main.cpp file are for error handling, MPI configuration and other special features.

In the constructor of the LAMMPS class instance the basic LAMMPS class hierachy
is created as shown in :ref:`class-topology`.  While processing the input further
class instances are created, or deleted, or replaced and specific member functions
of specific classes are called to trigger actions like creating atoms, computing
forces, computing properties, propagating the system, or writing output.


Inheritance and Compositing
===========================

Polymorphism
============


I/O and output formatting
^^^^^^^^^^^^^^^^^^^^^^^^^

Memory management
^^^^^^^^^^^^^^^^^


