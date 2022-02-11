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
if they are assigned to a pointer using the type of a base class.
Thanks to dynamic dispatch, LAMMPS can even use styles that are loaded
at runtime from a shared object file with the :doc:`plugin command <plugin>`.

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
or API of functions.

However, there are downsides to this. For example, calls virtual functions
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
   or default arguments for virtual functions.

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

Memory management
^^^^^^^^^^^^^^^^^
