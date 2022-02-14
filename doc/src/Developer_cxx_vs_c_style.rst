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

LAMMPS is designed to be an object oriented code, that is each
simulation is represented by an instance of the LAMMPS class.  When
running in parallel, of course, each MPI process will create such an
instance.  This can be seen in the ``main.cpp`` file where the core
steps of running a LAMMPS simulation are the following 3 lines of code:

.. code-block:: C++

    LAMMPS *lammps = new LAMMPS(argc, argv, lammps_comm);
    lammps->input->file();
    delete lammps;

The first line creates a LAMMPS class instance and passes the command
line arguments and the global communicator to its constructor.  The
second line tells the LAMMPS instance to process the input (either from
standard input or the provided input file) until the end.  And the third
line deletes that instance again.  The remainder of the main.cpp file
are for error handling, MPI configuration and other special features.

In the constructor of the LAMMPS class instance the basic LAMMPS class hierarchy
is created as shown in :ref:`class-topology`.  While processing the input further
class instances are created, or deleted, or replaced and specific member functions
of specific classes are called to trigger actions like creating atoms, computing
forces, computing properties, propagating the system, or writing output.

Compositing and Inheritance
===========================

LAMMPS makes extensive use of the object oriented programming (OOP)
principles of *compositing* and *inheritance*. Classes like the
``LAMMPS`` class are a **composite** containing pointers to instances of
other classes like ``Atom``, ``Comm``, ``Force``, ``Neighbor``,
``Modify``, and so on. Each of these classes implement certain
functionality by storing and manipulating data related to the simulation
and providing member functions that trigger certain actions.  Some of
those classes like ``Force`` and a composite again containing instances
of classes describing the force interactions or ``Modify`` containing
and calling fixes and computes. In most cases there is only one instance
of those member classes allowed, but in a few cases there can also be
multiple instances and the parent class is maintaining a list of the
pointers of instantiated classes.

Changing behavior or adjusting how LAMMPS handles the simulation is
implemented via **inheritance** where different variants of the
functionality are realized by creating *derived* classes that can share
common functionality in their base class and provide a consistent
interface where the derived classes replace (dummy or pure) functions in
the base class. The higher level classes can then call those methods of
the instantiated classes without having to know which specific derived
class variant was instantiated.  In the LAMMPS documentation those
derived classes are usually referred to a "styles", e.g.  pair styles,
fix styles, atom styles and so on.

This is the origin of the flexibility of LAMMPS and facilitates for
example to compute forces for very different non-bonded potential
functions by having different pair styles (implemented as different
classes derived from the ``Pair`` class) where the evaluation of the
potential function is confined to the implementation of the individual
classes.  Whenever a new :doc:`pair_style` or :doc:`bond_style` or
:doc:`comm_style` or similar command is processed in the LAMMPS input
any existing class instance is deleted and a new instance created in
it place.

Further code sharing is possible by creating derived classes from the
derived classes (for instance to implement an accelerated version of a
pair style) where then only a subset of the methods are replaced with
the accelerated versions.

Polymorphism
============

Polymorphism and dynamic dispatch are another OOP feature that play an
important part of how LAMMPS selects which code to execute.  In a nutshell,
this is a mechanism where the decision of which member function to call
from a class is determined at runtime and not when the code is compiled.
To enable it, the function has to be declared as ``virtual`` and all
corresponding functions in derived classes should be using the ``override``
property. Below is a brief example.

.. code-block:: c++

   class Base {
   public:
    virtual ~Base() = default;
    void call();
    void normal();
    virtual void poly();
   };

   void Base::call() {
    normal();
    poly();
   }

   class Derived : public Base {
   public:
    ~Derived() override = default;
    void normal();
    void poly() override;
   };

   // [....]

   Base *base1 = new Base();
   Base *base2 = new Derived();

   base1->call();
   base2->call();

The difference in behavior of the ``normal()`` and the ``poly()`` member
functions is in which of the two member functions is called when
executing `base1->call()` and `base2->call()`.  Without polymorphism, a
function within the base class will call only member functions within
the same scope, that is ``Base::call()`` will always call
``Base::normal()``.  But for the `base2->call()` the call for the
virtual member function will be dispatched to ``Derived::poly()``
instead.  This mechanism allows to always call functions within the
scope of the class type that was used to create the class instance, even
if they are assigned to a pointer using the type of a base class. This
is the desired behavior, and thanks to dynamic dispatch, LAMMPS can even
use styles that are loaded at runtime from a shared object file with the
:doc:`plugin command <plugin>`.

A special case of virtual functions are so-called pure functions. These
are virtual functions that are initialized to 0 in the class declaration
(see example below).

.. code-block:: c++

   class Base {
   public:
    virtual void pure() = 0;
   };

This has the effect that it will no longer be possible to create an instance
of the base class and that derived classes **must** implement these classes.
Many of the functions listed with the various styles in the section :doc:`Modify`
are such pure functions. The motivation for this is to define the interface
or API of functions but defer the implementation of those functionality to
the derived classes.

However, there are downsides to this. For example, calls to virtual functions
from within a constructor, will not be in the scope of the derived class and thus
it is good practice to either avoid calling them or to provide an explicit scope like
in ``Base::poly()``.  Furthermore, any destructors in classes containing
virtual functions should be declared virtual, too, so they are processed
in the expected order before types are removed from dynamic dispatch.

.. admonition:: Important Notes

   In order to be able to detect incompatibilities and to avoid unexpected
   behavior already at compile time, it is crucial that all member functions
   that are intended to replace a virtual or pure function use the ``override``
   property keyword.  For the same reason it should be avoided to use overloads
   or default arguments for virtual functions as they lead to confusion over
   which function is supposed to override which and which arguments need to be
   declared.

Style Factories
===============

In order to create class instances of the different styles, LAMMPS often
uses a programming pattern called `Factory`.  Those are functions that create
an instance of a specific derived class, say ``PairLJCut`` and return a pointer
to the type of the common base class of that style, ``Pair`` in this case.
To associate the factory function with the style keyword, an ``std::map``
class is used in which function pointers are indexed by their keyword
(for example "lj/cut" for ``PairLJCut`` and "morse" ``PairMorse``).
A couple of typedefs help to keep the code readable and a template function
is used to implement the actual factory functions for the individual classes.

I/O and output formatting
^^^^^^^^^^^^^^^^^^^^^^^^^

C-style stdio versus C++ style iostreams
========================================

LAMMPS chooses to use the "stdio" library of the standard C library for
reading from and writing to files and console instead of "iostreams" that were
introduced with C++.  This is mainly motivated by the better performance,
better control over formatting, and less effort to achieve specific formatting.

Since mixing "stdio" and "iostreams" can lead to unexpected behavior using
the latter is strongly discouraged.  Also output to the screen should not
use the predefined ``stdout`` FILE pointer, but rather the ``screen`` and
``logfile`` FILE pointers managed by the LAMMPS class.  Furthermore, output
should only be done by MPI rank 0 (``comm->me == 0``) and output that is
send to both ``screen`` and ``logfile`` should use the
:cpp:func:`utils::logmesg() convenience function <LAMMPS_NS::utils::logmesg>`.

Formatting with the {fmt} library
===================================

The LAMMPS source code includes a copy of the `{fmt} library
<https://fmt.dev>`_ which is preferred over formatting with the
"printf()" family of functions.  The primary reason is that it allows a
typesafe default format for any type of supported data.  This is
particularly useful for formatting integers of a given size (32-bit or
64-bit) which may require different format strings depending on compile
time settings or compilers/operating systems.  Furthermore, {fmt} gives
better performance, has more functionality, a familiar formatting syntax
that has similarities to ``format()`` in Python, and provides a facility
that can be used to integrate format strings and a variable number of
arguments into custom functions in a much simpler way that the varargs
mechanism of the C library.  Finally, {fmt} has been included into the
C++20 language standard, so changes to adopt it are future proof.

Formatted strings are frequently created by calling the
``fmt::format()`` function which will return a string as ``std::string``
class instance.  In contrast to the ``%`` placeholder in ``printf()``,
the {fmt} library uses ``{}`` to embed format descriptors.  In the
simplest case, no additional characters are needed as {fmt} will choose
the default format based on the data type of the argument. Alternatively
The ``fmt::print()`` function may be used instead of ``printf()`` or
``fprintf()``.  In addition, several LAMMPS output functions, that
originally accepted a single string as arguments have been overloaded to
accept a format string with optional arguments as well (e.g.
``Error::all()``, ``Error::one()``, ``utils::logmesg()``).

Summary of the {fmt} format syntax
==================================

The syntax of the format string is "{[<argument id>][:<format spec>]}",
where either the argument id or the format spec (separated by a colon
':') is optional.  The argument id is usually a number starting from 0
that is the index to the arguments following the format string.  By
default these are assigned in order (i.e. 0, 1, 2, 3, 4 etc.).  The most
common case for using argument id would be to use the same argument in
multiple places in the format string without having to provide it as an
argument multiple times. In LAMMPS the argument id is rarely used.

More common is the use of the format specifier, which starts with a
colon.  This may optionally be followed by a fill character (default is
' '). If provided, the fill character **must** be followed by an
alignment character ('<', '^', '>' for left, centered, or right
alignment (default)). The alignment character may be used without a fill
character. The next important format parameter would be the minimum
width, which may be followed by a dot '.'  and a precision for floating
point numbers. The final character in the format string would be an
indicator for the "presentation", i.e. 'd' for decimal presentation of
integers, 'x' for hexadecimal, 'o' for octal, 'c' for character
etc. This mostly follows the "printf()" scheme but without requiring an
additional length parameter to distinguish between different integer
widths. The {fmt} library will detect those and adapt the formatting
accordingly.  For floating point numbers there are correspondingly, 'g'
for generic presentation, 'e' for exponential presentation, and 'f' for
fixed point presentation.

Thus "{:8}" would represent *any* type argument using at least 8
characters; "{:<8}" would do this as left aligned, "{:^8}" as centered,
"{:>8}" as right aligned.  If a specific presentation is selected, the
argument type must be compatible or else the {fmt} formatting code will
throw an exception. Some format string examples are given below:

.. code-block:: C

   auto mesg = fmt::format("  CPU time: {:4d}:{:02d}:{:02d}\n", cpuh, cpum, cpus);
   mesg = fmt::format("{:<8s}| {:<10.5g} | {:<10.5g} | {:<10.5g} |{:6.1f} |{:6.2f}\n",
                      label, time_min, time, time_max, time_sq, tmp);
   utils::logmesg(lmp,"{:>6} = max # of 1-2 neighbors\n",maxall);
   utils::logmesg(lmp,"Lattice spacing in x,y,z = {:.8} {:.8} {:.8}\n",
                  xlattice,ylattice,zlattice);

A special feature of the {fmt} library is that format parameters like
the width or the precision may be also provided as arguments. In that
case a nested format is used where a pair of curly braces (with an
optional argument id) "{}" are used instead of the value, for example
"{:{}d}" will consume two integer arguments, the first will be the value
shown and the second the minimum width.

For more details and examples, please consult the `{fmt} syntax
documentation <https://fmt.dev/latest/syntax.html>`_ website.


Memory management
^^^^^^^^^^^^^^^^^
