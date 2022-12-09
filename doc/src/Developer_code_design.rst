Code design
-----------

This section explains some of the code design choices in LAMMPS with
the goal of helping developers write new code similar to the existing
code.  Please see the section on :doc:`Requirements for contributed
code <Modify_style>` for more specific recommendations and guidelines.
While that section is organized more in the form of a checklist for
code contributors, the focus here is on overall code design strategy,
choices made between possible alternatives, and discussing some
relevant C++ programming language constructs.

Historically, the basic design philosophy of the LAMMPS C++ code was a
"C with classes" style.  The motivation was to make it easy to modify
LAMMPS for people without significant training in C++ programming.
Data structures and code constructs were used that resemble the
previous implementation(s) in Fortran.  A contributing factor to this
choice also was that at the time, C++ compilers were often not mature
and some of the advanced features contained bugs or did not function
as the standard required.  There were also disagreements between
compiler vendors as to how to interpret the C++ standard documents.

However, C++ compilers have now advanced significantly.  In 2020 we
decided to to require the C++11 standard as the minimum C++ language
standard for LAMMPS.  Since then we have begun to also replace some of
the C-style constructs with equivalent C++ functionality, either from
the C++ standard library or as custom classes or functions, in order
to improve readability of the code and to increase code reuse through
abstraction of commonly used functionality.

.. note::

   Please note that as of spring 2022 there is still a sizable chunk
   of legacy code in LAMMPS that has not yet been refactored to
   reflect these style conventions in full.  LAMMPS has a large code
   base and many different contributors and there also is a hierarchy
   of precedence in which the code is adapted.  Highest priority has
   been the code in the ``src`` folder, followed by code in packages
   in order of their popularity and complexity (simpler code is
   adapted sooner), followed by code in the ``lib`` folder.  Source
   code that is downloaded from external packages or libraries during
   compilation is not subject to the conventions discussed here.

Object oriented code
^^^^^^^^^^^^^^^^^^^^

LAMMPS is designed to be an object oriented code.  Each simulation is
represented by an instance of the LAMMPS class.  When running in
parallel each MPI process creates such an instance.  This can be seen
in the ``main.cpp`` file where the core steps of running a LAMMPS
simulation are the following 3 lines of code:

.. code-block:: c++

    LAMMPS *lammps = new LAMMPS(argc, argv, lammps_comm);
    lammps->input->file();
    delete lammps;

The first line creates a LAMMPS class instance and passes the command
line arguments and the global communicator to its constructor.  The
second line triggers the LAMMPS instance to process the input (either
from standard input or a provided input file) until the simulation
ends.  The third line deletes the LAMMPS instance.  The remainder of
the main.cpp file has code for error handling, MPI configuration, and
other special features.

The basic LAMMPS class hierarchy which is created by the LAMMPS class
constructor is shown in :ref:`class-topology`.  When input commands
are processed, additional class instances are created, or deleted, or
replaced.  Likewise specific member functions of specific classes are
called to trigger actions such creating atoms, computing forces,
computing properties, time-propagating the system, or writing output.

Compositing and Inheritance
===========================

LAMMPS makes extensive use of the object oriented programming (OOP)
principles of *compositing* and *inheritance*. Classes like the
``LAMMPS`` class are a **composite** containing pointers to instances
of other classes like ``Atom``, ``Comm``, ``Force``, ``Neighbor``,
``Modify``, and so on.  Each of these classes implements certain
functionality by storing and manipulating data related to the
simulation and providing member functions that trigger certain
actions.  Some of those classes like ``Force`` are themselves
composites, containing instances of classes describing different force
interactions.  Similarly the ``Modify`` class contains a list of
``Fix`` and ``Compute`` classes.  If the input commands that
correspond to these classes include the word *style*, then LAMMPS
stores only a single instance of that class.  E.g. *atom_style*,
*comm_style*, *pair_style*, *bond_style*.  If the input command does
**not** include the word *style*, then there may be many instances of
that class defined, for example *region*, *fix*, *compute*, *dump*.

**Inheritance** enables creation of *derived* classes that can share
common functionality in their base class while providing a consistent
interface.  The derived classes replace (dummy or pure) functions in
the base class.  The higher level classes can then call those methods
of the instantiated classes without having to know which specific
derived class variant was instantiated.  In LAMMPS these derived
classes are often referred to as "styles", e.g.  pair styles, fix
styles, atom styles and so on.

This is the origin of the flexibility of LAMMPS.  For example pair
styles implement a variety of different non-bonded interatomic
potentials functions.  All details for the implementation of a
potential are stored and executed in a single class.

As mentioned above, there can be multiple instances of classes derived
from the ``Fix`` or ``Compute`` base classes.  They represent a
different facet of LAMMPS flexibility as they provide methods which
can be called at different points in time within a timestep, as
explained in `Developer_flow`.  This allows the input script to tailor
how a specific simulation is run, what diagnostic computations are
performed, and how the output of those computations is further
processed or output.

Additional code sharing is possible by creating derived classes from the
derived classes (e.g., to implement an accelerated version of a pair
style) where only a subset of the derived class methods are replaced
with accelerated versions.

Polymorphism
============

Polymorphism and dynamic dispatch are another OOP feature that play an
important role in how LAMMPS selects what code to execute.  In a
nutshell, this is a mechanism where the decision of which member
function to call from a class is determined at runtime and not when
the code is compiled.  To enable it, the function has to be declared
as ``virtual`` and all corresponding functions in derived classes
should use the ``override`` property. Below is a brief example.

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
functions is which of the two member functions is called when executing
`base1->call()` versus `base2->call()`.  Without polymorphism, a
function within the base class can only call member functions within the
same scope, that is ``Base::call()`` will always call
``Base::normal()``.  But for the `base2->call()` case the call of the
virtual member function will be dispatched to ``Derived::poly()``
instead.  This mechanism means that functions are called within the
scope of the class type that was used to *create* the class instance are
invoked; even if they are assigned to a pointer using the type of a base
class.  This is the desired behavior and this way LAMMPS can even use
styles that are loaded at runtime from a shared object file with the
:doc:`plugin command <plugin>`.

A special case of virtual functions are so-called pure functions.  These
are virtual functions that are initialized to 0 in the class declaration
(see example below).

.. code-block:: c++

   class Base {
   public:
    virtual void pure() = 0;
   };

This has the effect that an instance of the base class cannot be
created and that derived classes **must** implement these functions.
Many of the functions listed with the various class styles in the
section :doc:`Modify` are pure functions.  The motivation for this is
to define the interface or API of the functions but defer their
implementation to the derived classes.

However, there are downsides to this. For example, calls to virtual
functions from within a constructor, will not be in the scope of the
derived class and thus it is good practice to either avoid calling them
or to provide an explicit scope such as ``Base::poly()`` or
``Derived::poly()``.  Furthermore, any destructors in classes containing
virtual functions should be declared virtual too, so they will be
processed in the expected order before types are removed from dynamic
dispatch.

.. admonition:: Important Notes

   In order to be able to detect incompatibilities at compile time and
   to avoid unexpected behavior, it is crucial that all member functions
   that are intended to replace a virtual or pure function use the
   ``override`` property keyword.  For the same reason, the use of
   overloads or default arguments for virtual functions should be
   avoided as they lead to confusion over which function is supposed to
   override which and which arguments need to be declared.

Style Factories
===============

In order to create class instances for different styles, LAMMPS often
uses a programming pattern called `Factory`.  Those are functions that
create an instance of a specific derived class, say ``PairLJCut`` and
return a pointer to the type of the common base class of that style,
``Pair`` in this case.  To associate the factory function with the
style keyword, an ``std::map`` class is used with function pointers
indexed by their keyword (for example "lj/cut" for ``PairLJCut`` and
"morse" for ``PairMorse``).  A couple of typedefs help keep the code
readable and a template function is used to implement the actual
factory functions for the individual classes.  Below is an example
of such a factory function from the ``Force`` class as declared in
``force.h`` and implemented in ``force.cpp``.  The file ``style_pair.h``
is generated during compilation and includes all main header files
(i.e. those starting with ``pair_``) of pair styles and then the
macro ``PairStyle()`` will associate the style name "lj/cut"
with a factory function creating an instance of the ``PairLJCut``
class.

.. code-block:: c++

   // from force.h
   typedef Pair *(*PairCreator)(LAMMPS *);
   typedef std::map<std::string, PairCreator> PairCreatorMap;
   PairCreatorMap *pair_map;

   // from force.cpp
   template <typename S, typename T> static S *style_creator(LAMMPS *lmp)
   {
     return new T(lmp);
   }

   // [...]

   pair_map = new PairCreatorMap();

   #define PAIR_CLASS
   #define PairStyle(key, Class) (*pair_map)[#key] = &style_creator<Pair, Class>;
   #include "style_pair.h"
   #undef PairStyle
   #undef PAIR_CLASS

   // from pair_lj_cut.h

   #ifdef PAIR_CLASS
   PairStyle(lj/cut,PairLJCut);
   #else
   // [...]

Similar code constructs are present in other files like ``modify.cpp`` and
``modify.h`` or ``neighbor.cpp`` and ``neighbor.h``.  Those contain
similar macros and include ``style_*.h`` files for creating class instances
of styles they manage.


I/O and output formatting
^^^^^^^^^^^^^^^^^^^^^^^^^

C-style stdio versus C++ style iostreams
========================================

LAMMPS uses the "stdio" library of the standard C library for reading
from and writing to files and console instead of C++ "iostreams".
This is mainly motivated by better performance, better control over
formatting, and less effort to achieve specific formatting.

Since mixing "stdio" and "iostreams" can lead to unexpected
behavior. use of the latter is strongly discouraged.  Also output to
the screen should not use the predefined ``stdout`` FILE pointer, but
rather the ``screen`` and ``logfile`` FILE pointers managed by the
LAMMPS class.  Furthermore, output should generally only be done by
MPI rank 0 (``comm->me == 0``).  Output that is sent to both
``screen`` and ``logfile`` should use the :cpp:func:`utils::logmesg()
convenience function <LAMMPS_NS::utils::logmesg>`.

We also discourage the use of stringstreams because the bundled {fmt}
library and the customized tokenizer classes can provide the same
functionality in a cleaner way with better performance.  This also
helps maintain a consistent programming syntax with code from many
different contributors.

Formatting with the {fmt} library
===================================

The LAMMPS source code includes a copy of the `{fmt} library
<https://fmt.dev>`_ which is preferred over formatting with the
"printf()" family of functions.  The primary reason is that it allows
a typesafe default format for any type of supported data.  This is
particularly useful for formatting integers of a given size (32-bit or
64-bit) which may require different format strings depending on
compile time settings or compilers/operating systems.  Furthermore,
{fmt} gives better performance, has more functionality, a familiar
formatting syntax that has similarities to ``format()`` in Python, and
provides a facility that can be used to integrate format strings and a
variable number of arguments into custom functions in a much simpler
way than the varargs mechanism of the C library.  Finally, {fmt} has
been included into the C++20 language standard, so changes to adopt it
are future-proof.

Formatted strings are frequently created by calling the
``fmt::format()`` function which will return a string as a
``std::string`` class instance.  In contrast to the ``%`` placeholder
in ``printf()``, the {fmt} library uses ``{}`` to embed format
descriptors.  In the simplest case, no additional characters are
needed as {fmt} will choose the default format based on the data type
of the argument.  Otherwise the ``fmt::print()`` function may be
used instead of ``printf()`` or ``fprintf()``.  In addition, several
LAMMPS output functions, that originally accepted a single string as
argument have been overloaded to accept a format string with optional
arguments as well (e.g., ``Error::all()``, ``Error::one()``,
``utils::logmesg()``).

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

More common is the use of a format specifier, which starts with a colon.
This may optionally be followed by a fill character (default is ' '). If
provided, the fill character **must** be followed by an alignment
character ('<', '^', '>' for left, centered, or right alignment
(default)).  The alignment character may be used without a fill
character.  The next important format parameter would be the minimum
width, which may be followed by a dot '.' and a precision for floating
point numbers.  The final character in the format string would be an
indicator for the "presentation", i.e. 'd' for decimal presentation of
integers, 'x' for hexadecimal, 'o' for octal, 'c' for character etc.
This mostly follows the "printf()" scheme but without requiring an
additional length parameter to distinguish between different integer
widths.  The {fmt} library will detect those and adapt the formatting
accordingly.  For floating point numbers there are correspondingly, 'g'
for generic presentation, 'e' for exponential presentation, and 'f' for
fixed point presentation.

Thus "{:8}" would represent *any* type argument using at least 8
characters; "{:<8}" would do this as left aligned, "{:^8}" as centered,
"{:>8}" as right aligned.  If a specific presentation is selected, the
argument type must be compatible or else the {fmt} formatting code will
throw an exception. Some format string examples are given below:

.. code-block:: c++

   auto mesg = fmt::format("  CPU time: {:4d}:{:02d}:{:02d}\n", cpuh, cpum, cpus);
   mesg = fmt::format("{:<8s}| {:<10.5g} | {:<10.5g} | {:<10.5g} |{:6.1f} |{:6.2f}\n",
                      label, time_min, time, time_max, time_sq, tmp);
   utils::logmesg(lmp,"{:>6} = max # of 1-2 neighbors\n",maxall);
   utils::logmesg(lmp,"Lattice spacing in x,y,z = {:.8} {:.8} {:.8}\n",
                  xlattice,ylattice,zlattice);

which will create the following output lines:

.. parsed-literal::

     CPU time:    0:02:16
     Pair    | 2.0133     | 2.0133     | 2.0133     |   0.0 | 84.21
          4 = max # of 1-2 neighbors
     Lattice spacing in x,y,z = 1.6795962 1.6795962 1.6795962

Finally, a special feature of the {fmt} library is that format
parameters like the width or the precision may be also provided as
arguments. In that case a nested format is used where a pair of curly
braces (with an optional argument id) "{}" are used instead of the
value, for example "{:{}d}" will consume two integer arguments, the
first will be the value shown and the second the minimum width.

For more details and examples, please consult the `{fmt} syntax
documentation <https://fmt.dev/latest/syntax.html>`_ website.


Memory management
^^^^^^^^^^^^^^^^^

Dynamical allocation of small data and objects can be done with the
the C++ commands "new" and "delete/delete[].  Large data should use
the member functions of the ``Memory`` class, most commonly,
``Memory::create()``, ``Memory::grow()``, and ``Memory::destroy()``,
which provide variants for vectors, 2d arrays, 3d arrays, etc.
These can also be used for small data.

The use of ``malloc()``, ``calloc()``, ``realloc()`` and ``free()``
directly is strongly discouraged.  To simplify adapting legacy code
into the LAMMPS code base the member functions ``Memory::smalloc()``,
``Memory::srealloc()``, and ``Memory::sfree()`` are available, which
perform additional error checks for safety.

Use of these custom memory allocation functions is motivated by the
following considerations:

- memory allocation failures on *any* MPI rank during a parallel run
  will trigger an immediate abort of the entire parallel calculation
  instead of stalling it
- a failing "new" will trigger an exception which is also captured by
  LAMMPS and triggers a global abort
- allocation of multi-dimensional arrays will be done in a C compatible
  fashion but so that the storage of the actual data is stored in one
  large contiguous block.  Thus when MPI communication is needed,
  the data can be communicated directly (similar to Fortran arrays).
- the "destroy()" and "sfree()" functions may safely be called on NULL
  pointers
- the "destroy()" functions will nullify the pointer variables making
  "use after free" errors easy to detect
- it is possible to use a larger than default memory alignment (not on
  all operating systems, since the allocated storage pointers must be
  compatible with ``free()`` for technical reasons)

In the practical implementation of code this means that any pointer
variables that are class members should be initialized to a
``nullptr`` value in their respective constructors.  That way it is
safe to call ``Memory::destroy()`` or ``delete[]`` on them before
*any* allocation outside the constructor.  This helps prevent memory
leaks.
