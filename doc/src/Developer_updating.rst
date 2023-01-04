Notes for updating code written for older LAMMPS versions
---------------------------------------------------------

This section documents how C++ source files that are available *outside
of the LAMMPS source distribution* (e.g. in external USER packages or as
source files provided as a supplement to a publication) that are written
for an older version of LAMMPS and thus need to be updated to be
compatible with the current version of LAMMPS.  Due to the active
development of LAMMPS it is likely to always be incomplete.  Please
contact developers@lammps.org in case you run across an issue that is not
(yet) listed here.  Please also review the latest information about the
LAMMPS :doc:`programming style conventions <Modify_style>`, especially
if you are considering to submit the updated version for inclusion into
the LAMMPS distribution.

Available topics in mostly chronological order are:

- `Setting flags in the constructor`_
- `Rename of pack/unpack_comm() to pack/unpack_forward_comm()`_
- `Use ev_init() to initialize variables derived from eflag and vflag`_
- `Use utils::numeric() functions instead of force->numeric()`_
- `Use utils::open_potential() function to open potential files`_
- `Simplify customized error messages`_
- `Use of "override" instead of "virtual"`_
- `Simplified and more compact neighbor list requests`_
- `Split of fix STORE into fix STORE/GLOBAL and fix STORE/PERATOM`_
- `Use Output::get_dump_by_id() instead of Output::find_dump()`_
- `Refactored grid communication using Grid3d/Grid2d classes instead of GridComm`_

----

Setting flags in the constructor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As LAMMPS gains additional functionality, new flags may need to be set
in the constructor or a class to signal compatibility with such features.
Most of the time the defaults are chosen conservatively, but sometimes
the conservative choice is the uncommon choice, and then those settings
need to be made when updating code.

Pair styles:

  - ``manybody_flag``: set to 1 if your pair style is not pair-wise additive
  - ``restartinfo``: set to 0 if your pair style does not store data in restart files


Rename of pack/unpack_comm() to pack/unpack_forward_comm()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 8Aug2014

In this change set the functions to pack data into communication buffers
and to unpack data from communication buffers for :doc:`forward
communications <Developer_comm_ops>` were renamed from ``pack_comm()``
and ``unpack_comm()`` to ``pack_forward_comm()`` and
``unpack_forward_comm()``, respectively.  Also the meaning of the return
value of these functions was changed: rather than returning the number
of items per atom stored in the buffer, now the total number of items
added (or unpacked) needs to be returned.  Here is an example from the
`PairEAM` class.  Of course the member function declaration in corresponding
header file needs to be updated accordingly.

Old:

.. code-block:: c++

   int PairEAM::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
   {
     int m = 0;
     for (int i = 0; i < n; i++) {
       int j = list[i];
       buf[m++] = fp[j];
     }
     return 1;
   }

New:

.. code-block:: c++

   int PairEAM::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
   {
     int m = 0;
     for (int i = 0; i < n; i++) {
       int j = list[i];
       buf[m++] = fp[j];
     }
     return m;
   }

.. note::

   Because the various "pack" and "unpack" functions are defined in the
   respective base classes as dummy functions doing nothing, and because
   of the the name mismatch the custom versions in the derived class
   will no longer be called, there will be no compilation error when
   this change is not applied.  Only calculations will suddenly produce
   incorrect results because the required forward communication calls
   will cease to function correctly.

Use ev_init() to initialize variables derived from eflag and vflag
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 29Mar2019

There are several variables that need to be initialized based on
the values of the "eflag" and "vflag" variables and since sometimes
there are new bits added and new variables need to be set to 1 or 0.
To make this consistent, across all styles, there is now an inline
function ``ev_init(eflag, vflag)`` that makes those settings
consistently and calls either ``ev_setup()`` or ``ev_unset()``.
Example from a pair style:

Old:

.. code-block:: c++

   if (eflag || vflag) ev_setup(eflag, vflag);
   else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

New:

.. code-block:: c++

   ev_init(eflag, vflag);

Not applying this change will not cause a compilation error, but
can lead to inconsistent behavior and incorrect tallying of
energy or virial.

Use utils::numeric() functions instead of force->numeric()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 18Sep2020

The "numeric()" conversion functions (including "inumeric()",
"bnumeric()", and "tnumeric()") have been moved from the Force class to
the utils namespace.  Also they take an additional argument that selects
whether the ``Error::all()`` or ``Error::one()`` function should be
called in case of an error.  The former should be used when *all* MPI
processes call the conversion function and the latter *must* be used
when they are called from only one or a subset of the MPI processes.

Old:

.. code-block:: c++

    val = force->numeric(FLERR, arg[1]);
    num = force->inumeric(FLERR, arg[2]);

New:

.. code-block:: c++

    val = utils::numeric(FLERR, true, arg[1], lmp);
    num = utils::inumeric(FLERR, false, arg[2], lmp);

.. seealso::

   :cpp:func:`utils::numeric() <LAMMPS_NS::utils::numeric>`,
   :cpp:func:`utils::inumeric() <LAMMPS_NS::utils::inumeric>`,
   :cpp:func:`utils::bnumeric() <LAMMPS_NS::utils::bnumeric>`,
   :cpp:func:`utils::tnumeric() <LAMMPS_NS::utils::tnumeric>`

Use utils::open_potential() function to open potential files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 18Sep2020

The :cpp:func:`utils::open_potential()
<LAMMPS_NS::utils::open_potential>` function must be used to replace
calls to ``force->open_potential()`` and should be used to replace
``fopen()`` for opening potential files for reading.  The custom
function does three additional steps compared to ``fopen()``: 1) it will
try to parse the ``UNITS:`` and ``DATE:`` metadata will stop with an
error on a units mismatch and will print the date info, if present, in
the log file; 2) for pair styles that support it, it will set up
possible automatic unit conversions based on the embedded unit
information and LAMMPS' current units setting; 3) it will not only try
to open a potential file at the given path, but will also search in the
folders listed in the ``LAMMPS_POTENTIALS`` environment variable.  This
allows to keep potential files in a common location instead of having to
copy them around for simulations.

Old:

.. code-block:: c++

   fp = force->open_potential(filename);
   fp = fopen(filename, "r");

New:

.. code-block:: c++

   fp = utils::open_potential(filename, lmp);

Simplify customized error messages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 14May2021

Aided by features of the bundled {fmt} library, error messages now
can have a variable number of arguments and the string will be interpreted
as a {fmt} style format string so that custom error messages can be
easily customized without having to use temporary buffers and ``sprintf()``.
Example:

Old:

.. code-block:: c++

   if (fptr == NULL) {
     char str[128];
     sprintf(str,"Cannot open AEAM potential file %s",filename);
     error->one(FLERR,str);
   }

New:

.. code-block:: c++

   if (fptr == nullptr)
     error->one(FLERR, "Cannot open AEAM potential file {}: {}", filename, utils::getsyserror());

Use of "override" instead of "virtual"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 17Feb2022

Since LAMMPS requires C++11 we switched to use the "override" keyword
instead of "virtual" to indicate polymorphism in derived classes.  This
allows the C++ compiler to better detect inconsistencies when an
override is intended or not.  Please note that "override" has to be
added to **all** polymorph functions in derived classes and "virtual"
*only* to the function in the base class (or the destructor).  Here is
an example from the ``FixWallReflect`` class:

Old:

.. code-block:: c++

   FixWallReflect(class LAMMPS *, int, char **);
   virtual ~FixWallReflect();
   int setmask();
   void init();
   void post_integrate();

New:

.. code-block:: c++

   FixWallReflect(class LAMMPS *, int, char **);
   ~FixWallReflect() override;
   int setmask() override;
   void init() override;
   void post_integrate() override;

This change set will neither cause a compilation failure, nor will it
change functionality, but if you plan to submit the updated code for
inclusion into the LAMMPS distribution, it will be requested for achieve
a consistent :doc:`programming style <Modify_style>`.

Simplified function names for forward and reverse communication
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 24Mar2022

Rather then using the function name to distinguish between the different
forward and reverse communication functions for styles, LAMMPS now uses
the type of the "this" pointer argument.

Old:

.. code-block:: c++

   comm->forward_comm_pair(this);
   comm->forward_comm_fix(this);
   comm->forward_comm_compute(this);
   comm->forward_comm_dump(this);
   comm->reverse_comm_pair(this);
   comm->reverse_comm_fix(this);
   comm->reverse_comm_compute(this);
   comm->reverse_comm_dump(this);

New:

.. code-block:: c++

   comm->forward_comm(this);
   comm->reverse_comm(this);

This change is **required** or else the code will not compile.

Simplified and more compact neighbor list requests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 24Mar2022

This change set reduces the amount of code required to request a
neighbor list.  It enforces consistency and no longer requires to change
internal data of the request.  More information on neighbor list
requests can be :doc:`found here <Developer_notes>`. Example from the
``ComputeRDF`` class:

Old:

.. code-block:: c++

   int irequest = neighbor->request(this,instance_me);
   neighbor->requests[irequest]->pair = 0;
   neighbor->requests[irequest]->compute = 1;
   neighbor->requests[irequest]->occasional = 1;
   if (cutflag) {
     neighbor->requests[irequest]->cut = 1;
     neighbor->requests[irequest]->cutoff = mycutneigh;
   }

New:

.. code-block:: c++

   auto req = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
   if (cutflag) req->set_cutoff(mycutneigh);

Public access to the ``NeighRequest`` class data members has been
removed so this update is **required** to avoid compilation failure.

Split of fix STORE into fix STORE/GLOBAL and fix STORE/PERATOM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 15Sep2022

This change splits the GLOBAL and PERATOM modes of fix STORE into two
separate fixes STORE/GLOBAL and STORE/PERATOM.  There was very little
shared code between the two fix STORE modes and the two different code
paths had to be prefixed with if statements.  Furthermore, some flags
were used differently in the two modes leading to confusion.  Splitting
the code into two fix styles, makes it more easily maintainable.  Since
these are internal fixes, there is no user visible change.

Old:

.. code-block:: c++

   #include "fix_store.h"

   FixStore *fix = dynamic_cast<FixStore *>(
      modify->add_fix(fmt::format("{} {} STORE peratom 1 13",id_pole,group->names[0]));

   FixStore *fix = dynamic_cast<FixStore *>(modify->get_fix_by_id(id_pole));

New:

.. code-block:: c++

   #include "fix_store_peratom.h"

   FixStorePeratom *fix = dynamic_cast<FixStorePeratom *>(
      modify->add_fix(fmt::format("{} {} STORE/PERATOM 1 13",id_pole,group->names[0]));

   FixStorePeratom *fix = dynamic_cast<FixStorePeratom *>(modify->get_fix_by_id(id_pole));

Old:

.. code-block:: c++

   #include "fix_store.h"

   FixStore *fix = dynamic_cast<FixStore *>(
      modify->add_fix(fmt::format("{} {} STORE global 1 1",id_fix,group->names[igroup]));

   FixStore *fix = dynamic_cast<FixStore *>(modify->get_fix_by_id(id_fix));

New:

.. code-block:: c++

   #include "fix_store_global.h"

   FixStoreGlobal *fix = dynamic_cast<FixStoreGlobal *>(
      modify->add_fix(fmt::format("{} {} STORE/GLOBAL 1 1",id_fix,group->names[igroup]));

   FixStoreGlobal *fix = dynamic_cast<FixStoreGlobal *>(modify->get_fix_by_id(id_fix));

This change is **required** or else the code will not compile.

Use Output::get_dump_by_id() instead of Output::find_dump()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 15Sep2022

The accessor function to individual dump style instances has been changed
from ``Output::find_dump()`` returning the index of the dump instance in
the list of dumps to ``Output::get_dump_by_id()`` returning a pointer to
the dump directly.  Example:

Old:

.. code-block:: c++

   int idump = output->find_dump(arg[iarg+1]);
   if (idump < 0)
     error->all(FLERR,"Dump ID in hyper command does not exist");
   memory->grow(dumplist,ndump+1,"hyper:dumplist");
   dumplist[ndump++] = idump;

   [...]

   if (dumpflag)
     for (int idump = 0; idump < ndump; idump++)
       output->dump[dumplist[idump]]->write();

New:

.. code-block:: c++

   auto idump = output->get_dump_by_id(arg[iarg+1]);
   if (!idump) error->all(FLERR,"Dump ID {} in hyper command does not exist", arg[iarg+1]);
   dumplist.emplace_back(idump);

   [...]

   if (dumpflag) for (auto idump : dumplist) idump->write();

This change is **required** or else the code will not compile.

Refactored grid communication using Grid3d/Grid2d classes instead of GridComm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 22Dec2022

The ``GridComm`` class was for creating and communicating distributed
grids was replaced by the ``Grid3d`` class with added functionality.
A ``Grid2d`` class was also added for additional flexibility.

The new functionality and commands using the two grid classes are
discussed on the following documentation pages:

- :doc:`Howto_grid`
- :doc:`Developer_grid`

If you have custom LAMMPS code, which uses the GridComm class, here are some notes
on how to adapt it for using the Grid3d class.

(1) The constructor has changed to allow the ``Grid3d`` / ``Grid2d``
    classes to partition the global grid across processors, both for
    owned and ghost grid cells.  Previously any class which called
    ``GridComm`` performed the partitioning itself and that information
    was passed in the ``GridComm::GridComm()`` constructor.  There are
    several "set" functions which can be called to alter how ``Grid3d``
    / ``Grid2d`` perform the partitioning.  They should be sufficient
    for most use cases of the grid classes.

(2) The partitioning is triggered by the ``setup_grid()`` method.

(3) The ``setup()`` method of the ``GridComm`` class has been replaced
    by the ``setup_comm()`` method in the new grid classes.  The syntax
    for the ``forward_comm()`` and ``reverse_comm()`` methods is
    slightly altered as is the syntax of the associated pack/unpack
    callback methods.  But the functionality of these operations is the
    same as before.

(4) The new ``Grid3d`` / ``Grid2d`` classes have additional
    functionality for dynamic load-balancing of grids and their
    associated data across processors.  This did not exist in the
    ``GridComm`` class.

This and more is explained in detail on the :doc:`Developer_grid` page.
The following LAMMPS source files can be used as illustrative examples
for how the new grid classes are used by computes, fixes, and various
KSpace solvers which use distributed FFT grids:

- ``src/fix_ave_grid.cpp``
- ``src/compute_property_grid.cpp``
- ``src/EXTRA-FIX/fix_ttm_grid.cpp``
- ``src/KSPACE/pppm.cpp``

This change is **required** or else the code will not compile.
